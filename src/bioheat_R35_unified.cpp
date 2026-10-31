// MFEM 4.10 / MPI-HYPRE: UNIFIED R35 BIOHEAT + ARRHENIUS SOLVER.
//
// ONE source combines both previously tested workflows:
//   (1) Mentor/reference Pennes + Arrhenius (10 W by default)
//   (2) Power/duration sensitivity (--power P, --tend t), with adjustable
//       blood perfusion, kinetic coefficients and modelled damage threshold.
//
// The validated 10 W DG0 Q source is matched to 106/107 tetrahedra by
// geometric centroid and volume. Original Q power is checked BEFORE any
// scaling. EM is not re-solved: Q(P)=P/10 * Q(10 W) assumes linear EM and
// temperature-independent properties. Each runtime writes energy and damage
// CSV plus time-dependent ParaView PVD. No explicit blood vessel geometry.
// All external thermal boundaries are adiabatic. Neither temperature nor
// Arrhenius thresholds demonstrate clinically validated necrosis.
//
// BACKWARD COMPATIBILITY: output basename is unchanged from the previously
// validated Arrhenius and power-sweep executables. Run tests in a separate
// working directory to avoid overwriting any reference CSV/PVD output.

#include "mfem.hpp"
#include <mpi.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef MFEM_USE_MPI
#error "This program requires MFEM_USE_MPI=ON"
#endif

namespace
{
constexpr int LIVER = 106, TUMOR = 107;
constexpr double DEFAULT_EXPECTED_POWER = 9.0049672684; // W, validated case
constexpr double HASH_GRID = 1.0e-6; // metre = 1 micrometre grid
constexpr double CENTROID_TOL = 2.0e-7; // 0.2 micrometre (must be unique)
constexpr double VOLUME_REL_TOL = 2.0e-3; // 0.2%
constexpr double GAS_R = 8.31446261815324; // J/(mol K)
constexpr double MAX_OMEGA = 1.0e12; // cap to avoid overflow; theta is then effectively 1
constexpr const char *OUT = "bioheat_R35_10W_arrhenius";


struct Input
{
   std::string mesh, qcsv;
   double tend=300.0, dt=1.0, omega106=0.0036, omega107=0.0036;
   double expectedQ=DEFAULT_EXPECTED_POWER; // validation of ORIGINAL 10 W input CSV
   double powerW=10.0; // nominal forward/incident EM power in W; input Q is for 10 W
   double initC=37.0, bloodC=37.0;
   int output_every=15;
   // Pilot kinetic parameters (literature for liver + liver tumor, not patient-calibrated).
   // NOTE: Some secondary sources misprint A_tumor exponent; check primary reference.
   double A106=7.39e39, Ea106=2.577e5; // 1/s, J/mol
   double A107=3.247e43, Ea107=2.814e5;
   double damage_threshold=0.98; // theta_d >= 0.98 <=> Omega >= -ln(0.02)

};

void Usage()
{
   std::cerr << "Usage: mpirun -np 2 ./bioheat_R35_unified MESH.msh Q_Pinc10W_cells.csv"
                " [--tend 300] [--dt 1] [--omega106 0.0036] [--omega107 0.0036]"
                " [--save-every 15] [--expected-q 9.0049672684] [--power 10]\n"
                "  [--A106 7.39e39] [--Ea106 2.577e5] [--A107 3.247e43] [--Ea107 2.814e5]\n"
                "  [--damage-threshold 0.98] (modelled fraction, not clinical diagnosis)\n"
                "For no perfusion: --omega106 0 --omega107 0\n";
}

// Use an explicit power label in all NEW output names; keep original 10 W
// filenames for compatibility with already-verified sensitivity scripts.
std::string PowerLabel(double powerW)
{
   std::ostringstream ss;
   ss << std::fixed << std::setprecision(3) << powerW;
   std::string result=ss.str();
   while(result.find('.')!=std::string::npos && !result.empty() && result.back()=='0')
      result.pop_back();
   if(!result.empty() && result.back()=='.') result.pop_back();
   std::replace(result.begin(),result.end(),'.','p');
   return result;
}

bool Parse(int argc, char **argv, Input &in)
{
   if(argc<3 || (argc-3)%2!=0) {return false;}
   in.mesh=argv[1]; in.qcsv=argv[2];
   try
   {
      for(int i=3;i<argc;i+=2)
      {
         const std::string k=argv[i], s=argv[i+1];
         if(k=="--tend") in.tend=std::stod(s);
         else if(k=="--dt") in.dt=std::stod(s);
         else if(k=="--omega106") in.omega106=std::stod(s);
         else if(k=="--omega107") in.omega107=std::stod(s);
         else if(k=="--save-every") in.output_every=std::stoi(s);
         else if(k=="--expected-q") in.expectedQ=std::stod(s);
         else if(k=="--power") in.powerW=std::stod(s);
         else if(k=="--A106") in.A106=std::stod(s);
         else if(k=="--Ea106") in.Ea106=std::stod(s);
         else if(k=="--A107") in.A107=std::stod(s);
         else if(k=="--Ea107") in.Ea107=std::stod(s);
         else if(k=="--damage-threshold") in.damage_threshold=std::stod(s);
         else return false;
      }
   }
   catch (const std::exception&) {return false;}
   return std::isfinite(in.tend) && in.tend>0 && std::isfinite(in.dt) && in.dt>0 &&
          std::isfinite(in.omega106) && in.omega106>=0 &&
          std::isfinite(in.omega107) && in.omega107>=0 &&
          std::isfinite(in.expectedQ) && in.expectedQ>0 &&
          std::isfinite(in.powerW) && in.powerW>0 && in.powerW<=1000.0 &&
          in.output_every>=1 &&
          std::isfinite(in.A106) && in.A106>0 && std::isfinite(in.Ea106) && in.Ea106>0 &&
          std::isfinite(in.A107) && in.A107>0 && std::isfinite(in.Ea107) && in.Ea107>0 &&
          std::isfinite(in.damage_threshold) && in.damage_threshold>0 && in.damage_threshold<1 &&
          std::abs(std::round(in.tend/in.dt)*in.dt-in.tend) < 1e-9*in.tend;
}

struct Key
{
   std::int64_t x,y,z;
   bool operator==(const Key &rhs) const {return x==rhs.x&&y==rhs.y&&z==rhs.z;}
};
struct KeyHash
{
   size_t operator()(const Key &k) const
   {
      size_t h=1469598103934665603ull;
      for (std::int64_t v : {k.x,k.y,k.z})
      {
         const size_t a=std::hash<std::int64_t>{}(v);
         h^=a; h*=1099511628211ull;
      }
      return h;
   }
};
Key MakeKey(const std::array<double,3> &p)
{
   return {static_cast<std::int64_t>(std::llround(p[0]/HASH_GRID)),
           static_cast<std::int64_t>(std::llround(p[1]/HASH_GRID)),
           static_cast<std::int64_t>(std::llround(p[2]/HASH_GRID))};
}
struct Record
{
   std::array<double,3> c;
   double volume,q;
};
using QTable = std::unordered_map<Key,std::vector<Record>,KeyHash>;

bool ReadQ(const std::string &filename,QTable &table, size_t &count)
{
   std::ifstream file(filename);
   if(!file) {std::cerr<<"Cannot read Q CSV: "<<filename<<'\n';return false;}
   std::string line;
   if(!std::getline(file,line) || line.find("cx_m")!=0 || line.find("q_W_m3")==std::string::npos)
   {
      std::cerr<<"Invalid Q CSV header; use export_q_from_paraview.py\n";
      return false;
   }
   count=0;
   while(std::getline(file,line))
   {
      if(line.empty())continue;
      std::stringstream ss(line);
      std::string token;
      double v[5];
      for(int j=0;j<5;++j)
      {
         if(!std::getline(ss,token,','))return false;
         try {v[j]=std::stod(token);} catch(const std::exception&) {return false;}
         if(!std::isfinite(v[j]))return false;
      }
      if(v[3]<=0 || v[4]<0)return false;
      Record rec={{v[0],v[1],v[2]},v[3],v[4]};
      table[MakeKey(rec.c)].push_back(rec);
      ++count;
   }
   return count>0;
}

void ElementCenterAndVolume(mfem::Mesh &mesh, int e,
                            std::array<double,3> &center, double &vol)
{
   mfem::Array<int> ids;
   mesh.GetElementVertices(e,ids);
   if(ids.Size()!=4)throw std::runtime_error("Non-tetrahedral thermal cell");
   center={0.,0.,0.};
   for(int j=0;j<4;j++)
   {
      const auto *v=mesh.GetVertex(ids[j]);
      for(int k=0;k<3;k++)center[k]+=0.25*v[k];
   }
   auto *tr=mesh.GetElementTransformation(e);
   const auto &ip=mfem::Geometries.GetCenter(mfem::Geometry::TETRAHEDRON);
   tr->SetIntPoint(&ip);
   vol=tr->Weight()/6.0;
   if(!(vol>0.0)) throw std::runtime_error("Negative/degenerate thermal element volume");
}

const Record *FindQ(const QTable &table, const std::array<double,3> &c, double vol)
{
   const Key key=MakeKey(c);
   const Record* selected=nullptr;
   for(std::int64_t dx=-1;dx<=1;dx++)
   for(std::int64_t dy=-1;dy<=1;dy++)
   for(std::int64_t dz=-1;dz<=1;dz++)
   {
      const auto it=table.find({key.x+dx,key.y+dy,key.z+dz});
      if(it==table.end())continue;
      for(const auto &r: it->second)
      {
         double d2=0;
         for(int k=0;k<3;k++)d2+=(r.c[k]-c[k])*(r.c[k]-c[k]);
         if(d2>=CENTROID_TOL*CENTROID_TOL)continue;
         if(std::abs(r.volume-vol)>VOLUME_REL_TOL*vol)continue;
         if(selected)throw std::runtime_error("Ambiguous Q match; abort rather than mis-map cells");
         selected=&r;
      }
   }
   return selected;
}

struct Sum
{
   double q106=0,q107=0,volume106=0,volume107=0;
   int n106=0,n107=0;
};

Sum LoadSource(mfem::ParMesh &mesh, mfem::ParFiniteElementSpace &q_fes,
               mfem::ParGridFunction &q, const QTable &table)
{
   mfem::Array<int> dofs;
   Sum s;
   q=0.0;
   for(int i=0;i<mesh.GetNE();i++)
   {
      const int attr=mesh.GetAttribute(i);
      if(attr!=LIVER && attr!=TUMOR)
         throw std::runtime_error("Unexpected thermal material attribute");
      std::array<double,3> c;
      double vol;
      ElementCenterAndVolume(mesh,i,c,vol);
      const Record *r=FindQ(table,c,vol);
      if(!r)throw std::runtime_error("Missing heat source for a tissue tetrahedron");
      q_fes.GetElementDofs(i,dofs);
      if(dofs.Size()!=1 || dofs[0]<0)throw std::runtime_error("DG0 mapping is not 1 DOF per tetrahedron");
      q(dofs[0])=r->q;
      if(attr==LIVER){s.q106+=r->q*vol;s.volume106+=vol;s.n106++;}
      else{s.q107+=r->q*vol;s.volume107+=vol;s.n107++;}
   }
   return s;
}

void AllreduceSum(MPI_Comm comm, Sum &s)
{
   double v[4]={s.q106,s.q107,s.volume106,s.volume107};
   double w[4]={};
   MPI_Allreduce(v,w,4,MPI_DOUBLE,MPI_SUM,comm);
   s.q106=w[0];s.q107=w[1];s.volume106=w[2];s.volume107=w[3];
   int n[2]={s.n106,s.n107}, counts[2];
   MPI_Allreduce(n,counts,2,MPI_INT,MPI_SUM,comm);
   s.n106=counts[0];s.n107=counts[1];
}

// Damage is integrated in time independently at tetrahedral quadrature points.
// Use trapezoid rule in TIME and a positive 2nd-order tetrahedral rule in SPACE.
// The spatial integration point values are retained between thermal time steps.
// This avoids ambiguous material assignment of H1 nodes on the 106/107 interface.
// DG0 outputs are volume-weighted element averages, not nodal damage values.
struct DamageStats
{
   double V106=0., V107=0.;            // volume classified damaged [m3]
   double thetaInt106=0., thetaInt107=0.; // int theta dV [m3]
   double maxOmega106=0., maxOmega107=0.;
};

double ArrheniusRate(double TC, double A, double Ea)
{
   const double TK=TC+273.15;
   if(!(TK>0) || !std::isfinite(TK))
      throw std::runtime_error("Invalid temperature for Arrhenius model");
   const double exponent=std::log(A)-Ea/(GAS_R*TK);
   // Once the rate becomes enormous, Omega is capped (fraction -> 1).
   if(exponent>=std::log(MAX_OMEGA))return MAX_OMEGA;
   if(exponent < -745.)return 0.0;
   return std::exp(exponent);
}

DamageStats AdvanceDamage(mfem::ParMesh &mesh,
                          const mfem::IntegrationRule &ir,
                          const mfem::ParGridFunction &rise_prev,
                          const mfem::ParGridFunction &rise_now,
                          double initialC, double dt,
                          double A106, double Ea106,
                          double A107, double Ea107,
                          double omega_cut,
                          std::vector<double> &omega_at_qp,
                          mfem::ParFiniteElementSpace &dgspace,
                          mfem::ParGridFunction &dgOmega,
                          mfem::ParGridFunction &dgTheta,
                          mfem::ParGridFunction &dgCoverage)
{
   DamageStats stats;
   mfem::Array<int> dofs;
   mfem::Vector prevVals,nowVals;
   const int nq=ir.GetNPoints();
   if(omega_at_qp.size()!=static_cast<size_t>(mesh.GetNE())*nq)
      throw std::runtime_error("Arrhenius quadrature storage mismatch");
   for(int e=0;e<mesh.GetNE();e++)
   {
      const int attr=mesh.GetAttribute(e);
      const double A=(attr==LIVER)?A106:A107;
      const double Ea=(attr==LIVER)?Ea106:Ea107;
      if(attr!=LIVER && attr!=TUMOR)
         throw std::runtime_error("Unexpected material in damage integration");
      // MFEM 4.10 GetValues is faster and avoids changing transformation
      // state repeatedly inside the quadrature loop.
      rise_prev.GetValues(e,ir,prevVals);
      rise_now.GetValues(e,ir,nowVals);
      if(prevVals.Size()!=nq || nowVals.Size()!=nq)
         throw std::runtime_error("Unexpected thermal quadrature field size");
      auto *tr=mesh.GetElementTransformation(e);
      double cellV=0, omegaV=0, thetaV=0,coveredV=0;
      for(int j=0;j<nq;j++)
      {
         const mfem::IntegrationPoint &ip=ir.IntPoint(j);
         tr->SetIntPoint(&ip);
         const double dV=ip.weight*tr->Weight();
         if(!(dV>0.) || !std::isfinite(dV))
            throw std::runtime_error("Invalid thermal quadrature weight");
         const double Tprev=initialC+prevVals(j);
         const double Tnow =initialC+nowVals(j);
         const double r0=ArrheniusRate(Tprev,A,Ea);
         const double r1=ArrheniusRate(Tnow,A,Ea);
         const size_t idx=static_cast<size_t>(e)*nq+j;
         const double increment=0.5*dt*(r0+r1);
         omega_at_qp[idx]=std::min(MAX_OMEGA,omega_at_qp[idx]+increment);
         const double om=omega_at_qp[idx];
         const double theta=-std::expm1(-om); // stable for small Omega
         cellV+=dV;
         omegaV+=dV*om;
         thetaV+=dV*theta;
         if(om>=omega_cut)coveredV+=dV;
         if(attr==LIVER)
         {
            stats.maxOmega106=std::max(stats.maxOmega106,om);
            stats.thetaInt106+=dV*theta;
            if(om>=omega_cut)stats.V106+=dV;
         }
         else
         {
            stats.maxOmega107=std::max(stats.maxOmega107,om);
            stats.thetaInt107+=dV*theta;
            if(om>=omega_cut)stats.V107+=dV;
         }
      }
      dgspace.GetElementDofs(e,dofs);
      if(dofs.Size()!=1 || dofs[0]<0 || !(cellV>0))
         throw std::runtime_error("DG0 output mapping is invalid");
      dgOmega(dofs[0])=omegaV/cellV;
      dgTheta(dofs[0])=thetaV/cellV;
      dgCoverage(dofs[0])=coveredV/cellV;
   }
   return stats;
}

void ReduceDamage(MPI_Comm comm, DamageStats &d)
{
   double local[4]={d.V106,d.V107,d.thetaInt106,d.thetaInt107};
   double global[4]={};
   MPI_Allreduce(local,global,4,MPI_DOUBLE,MPI_SUM,comm);
   d.V106=global[0];d.V107=global[1];
   d.thetaInt106=global[2];d.thetaInt107=global[3];
   double locmax[2]={d.maxOmega106,d.maxOmega107};
   double globmax[2]={};
   MPI_Allreduce(locmax,globmax,2,MPI_DOUBLE,MPI_MAX,comm);
   d.maxOmega106=globmax[0];d.maxOmega107=globmax[1];
}

double GlobalDot(MPI_Comm comm, const mfem::Vector &a, const mfem::Vector &b)
{
   double loc=a*b, result=0;
   MPI_Allreduce(&loc,&result,1,MPI_DOUBLE,MPI_SUM,comm);
   return result;
}

double GlobalMax(MPI_Comm comm, double x)
{
   double y;
   MPI_Allreduce(&x,&y,1,MPI_DOUBLE,MPI_MAX,comm);
   return y;
}

double GlobalMin(MPI_Comm comm, double x)
{
   double y;
   MPI_Allreduce(&x,&y,1,MPI_DOUBLE,MPI_MIN,comm);
   return y;
}
} // namespace

int main(int argc, char **argv)
{
   MPI_Init(&argc,&argv);
   const MPI_Comm comm=MPI_COMM_WORLD;
   int rank=0;
   MPI_Comm_rank(comm,&rank);
   Input in;
   if(!Parse(argc,argv,in))
   {
      if(rank==0)Usage();
      MPI_Finalize();
      return 1;
   }
   try
   {
      // Every rank loads the same source and mesh; no partition ordering reliance.
      QTable table;
      size_t csv_count=0;
      if(!ReadQ(in.qcsv,table,csv_count))throw std::runtime_error("Failed to load source CSV");
      std::ifstream inp(in.mesh,std::ios::binary);
      if(!inp)throw std::runtime_error("Cannot read input .msh mesh");
      mfem::Mesh full(inp,1,0);
      if(full.Dimension()!=3 || full.SpaceDimension()!=3)
         throw std::runtime_error("Expected 3D tetrahedral Gmsh mesh in metres");
      mfem::Array<int> tissue_attrs;
      tissue_attrs.Append(LIVER);tissue_attrs.Append(TUMOR);
      mfem::SubMesh serial_tissue=mfem::SubMesh::CreateFromDomain(full,tissue_attrs);
      mfem::ParMesh mesh(comm,serial_tissue);
      if(rank==0)
      {
         std::cout << "=== MFEM 4.10 PARALLEL BIOHEAT (TISSUE 106/107) ===\n"
                   << "Input mesh: " << in.mesh << "\nQ CSV: " << in.qcsv << " ("
                   << csv_count << " EM tetrahedra)\n";
      }

      mfem::L2_FECollection qfec(0,3);
      mfem::ParFiniteElementSpace qfes(&mesh,&qfec);
      mfem::ParGridFunction q(&qfes);
      Sum sums=LoadSource(mesh,qfes,q,table);
      // CSV hash table no longer needed; release its memory before solving.
      table.clear(); table.rehash(0);
      AllreduceSum(comm,sums);
      const double pref=sums.q106+sums.q107;
      // Always check the ORIGINAL 10 W CSV BEFORE scaling it.
      if(!std::isfinite(pref) || pref<=0 ||
         std::abs(pref-in.expectedQ)>1e-3*in.expectedQ)
         throw std::runtime_error("Original 10 W Q mismatches validated reference by >0.1%; stop");
      const double scale=in.powerW/10.0;
      q *= scale; // scale all 106/107 DG0 cell power densities [W/m3]
      sums.q106 *= scale;
      sums.q107 *= scale;
      const double Psrc=sums.q106+sums.q107;
      const bool original10=(std::abs(in.powerW-10.0)<1.0e-12);
      const std::string powerLabel=PowerLabel(in.powerW);
      const std::string outName=original10 ? std::string(OUT)
          : ("bioheat_R35_Pinc"+powerLabel+"W_arrhenius");
      const std::string qFieldName=original10 ? std::string("Q_Pinc_10W_W_m3")
          : ("Q_Pinc_"+powerLabel+"W_W_m3");
      if(rank==0)
      {
         std::cout<<std::setprecision(12)
                  << "Thermal tetrahedra: 106=" << sums.n106
                  << " 107=" << sums.n107 << "\n"
                  << "Volume 106 [mm3]=" << 1e9*sums.volume106
                  << " 107 [mm3]=" << 1e9*sums.volume107 << "\n"
                  << "Original source integral at Pinc=10 W [W]="<<pref<<"\n"
                  << "Pinc [W]="<<in.powerW<<" ; Q scale="<<scale<<"\n"
                  << "Scaled source power 106 [W]="<<sums.q106<<" 107 [W]="<<sums.q107
                  <<" TOTAL [W]="<<Psrc<<'\n';
      }
      // Literature-based, temperature-independent pilot data.
      // Only tissue submesh participates. All thermal boundaries are insulated.
      // Physical continuity across 106/107 follows from continuous H1 trial space.
      mfem::Vector vc(107), vk(107), vp(107);
      vc=0.0;vk=0.0;vp=0.0;
      vc(LIVER-1)=1079.0*3540.0; // rho*Cp in J/(m^3 K)
      vc(TUMOR-1)=1040.0*3960.0;
      vk(LIVER-1)=0.52; // W/(m K)
      vk(TUMOR-1)=0.57;
      vp(LIVER-1)=1060.0*3639.0*in.omega106; // W/(m3 K)
      vp(TUMOR-1)=1060.0*3639.0*in.omega107;
      mfem::PWConstCoefficient heat_cap(vc), thermal_k(vk), perf(vp);
      mfem::GridFunctionCoefficient qcoeff(&q);
      mfem::H1_FECollection tfec(1,3);
      mfem::ParFiniteElementSpace fes(&mesh,&tfec);
      if(rank==0)std::cout<<"Temperature global true DOF: "<<fes.GlobalTrueVSize()<<'\n';
      mfem::ParBilinearForm mass(&fes);
      mass.AddDomainIntegrator(new mfem::MassIntegrator(heat_cap));
      mass.Assemble();mass.Finalize();
      mfem::ParBilinearForm stiffness(&fes);
      stiffness.AddDomainIntegrator(new mfem::DiffusionIntegrator(thermal_k));
      stiffness.AddDomainIntegrator(new mfem::MassIntegrator(perf));
      stiffness.Assemble();stiffness.Finalize();
      std::unique_ptr<mfem::HypreParMatrix> M(mass.ParallelAssemble());
      std::unique_ptr<mfem::HypreParMatrix> K(stiffness.ParallelAssemble());
      std::unique_ptr<mfem::HypreParMatrix> S(mfem::Add(1.0,*M,in.dt,*K));

      // Assemble source and lump-independent energy-perfusion integrals.
      mfem::ParLinearForm src(&fes), eform(&fes), pform(&fes);
      src.AddDomainIntegrator(new mfem::DomainLFIntegrator(qcoeff));
      eform.AddDomainIntegrator(new mfem::DomainLFIntegrator(heat_cap));
      pform.AddDomainIntegrator(new mfem::DomainLFIntegrator(perf));
      src.Assemble();eform.Assemble();pform.Assemble();
      std::unique_ptr<mfem::HypreParVector> F(src.ParallelAssemble());

      // Diagnostic: verify that every degree of freedom has a strictly
      // positive, finite diagonal in the implicit heat-system matrix.
      mfem::Vector sdiag;
      S->GetDiag(sdiag);
      double dmin_local = std::numeric_limits<double>::infinity();
      double dmax_local = 0.0;
      int bad_diag_local = 0;
      for (int i = 0; i < sdiag.Size(); ++i)
      {
         const double d = sdiag(i);
         if (!std::isfinite(d) || d <= 0.0) { ++bad_diag_local; }
         else
         {
            dmin_local = std::min(dmin_local, d);
            dmax_local = std::max(dmax_local, d);
         }
      }
      double dmin = 0.0, dmax = 0.0;
      int bad_diag = 0;
      MPI_Allreduce(&dmin_local, &dmin, 1, MPI_DOUBLE, MPI_MIN, comm);
      MPI_Allreduce(&dmax_local, &dmax, 1, MPI_DOUBLE, MPI_MAX, comm);
      MPI_Allreduce(&bad_diag_local, &bad_diag, 1, MPI_INT, MPI_SUM, comm);
      if (rank == 0)
      {
         std::cout << "Heat system S=M+dt(K+P): min(diag)=" << dmin
                   << ", max(diag)=" << dmax
                   << ", nonpositive/nonfinite diagonals=" << bad_diag
                   << '\n' << std::flush;
      }
      if (bad_diag != 0)
         throw std::runtime_error("Invalid heat-system matrix diagonal: stop before CG");

      // First isolate the BoomerAMG setup error using a simpler, SPD-compatible
      // Jacobi preconditioner. Do not assume this proves the matrix is SPD.
      mfem::HypreSmoother jacobi(*S, mfem::HypreSmoother::Jacobi);
      mfem::CGSolver cg(comm);
      cg.SetPreconditioner(jacobi);
      cg.SetOperator(*S);
      cg.SetRelTol(1e-11);
      cg.SetAbsTol(1e-13);
      cg.SetMaxIter(3000);
      cg.SetPrintLevel(0);
      if (rank == 0)
         std::cout << "Preconditioner: Jacobi (CG/HYPRE; existing validated configuration)\n"
                   << std::flush;

      mfem::ParGridFunction rise(&fes),temperature(&fes),rise_prev(&fes);
      rise=0.0;rise_prev=0.0;temperature=in.initC;
      mfem::Vector x(fes.GetTrueVSize()),rhs(fes.GetTrueVSize());
      x=0.0;
      const mfem::IntegrationRule &damageIR=mfem::IntRules.Get(mfem::Geometry::TETRAHEDRON,2);
      if(damageIR.GetNPoints()<4)
         throw std::runtime_error("Need positive multi-point tetrahedral integration for damage");
      const double omega_cut=-std::log1p(-in.damage_threshold);
      std::vector<double> omega_at_qp(static_cast<size_t>(mesh.GetNE())*damageIR.GetNPoints(),0.0);
      mfem::ParGridFunction dgOmega(&qfes), dgTheta(&qfes), dgCoverage(&qfes);
      dgOmega=0.0;dgTheta=0.0;dgCoverage=0.0;
      mfem::ParaViewDataCollection pv(outName,&mesh);
      pv.SetPrefixPath("rezultati");
      pv.SetLevelsOfDetail(1);
      pv.SetDataFormat(mfem::VTKFormat::BINARY);
      pv.RegisterField("Temperature_C",&temperature);
      pv.RegisterField("TemperatureRise_K",&rise);
      pv.RegisterField(qFieldName,&q);
      pv.RegisterField("Damage_Omega_DG0_avg",&dgOmega);
      pv.RegisterField("NecroticFraction_DG0_avg",&dgTheta);
      // Its integral over material 107 / tumor volume equals the CSV coverage.
      pv.RegisterField("DamageThresholdCoverage_DG0",&dgCoverage);
      auto save=[&](int step,double t)
      {
         pv.SetCycle(step);
         pv.SetTime(t);
         pv.Save();
         if(pv.Error()!=0)throw std::runtime_error("ParaView write failed");
      };
      std::ofstream csv,damage_csv;
      if(rank==0)
      {
         csv.open(outName+"_energy.csv");
         if(!csv)throw std::runtime_error("Unable to open CSV output");
         damage_csv.open(outName+"_damage.csv");
         if(!damage_csv)throw std::runtime_error("Unable to open damage CSV output");
         damage_csv<<"time_s,theta_threshold,Omega_threshold,damaged_106_cm3,damaged_107_cm3,"
                     "tumor_coverage_pct,healthy_damage_pct_of_106,healthy_damaged_cm3,"
                     "mean_theta106,mean_theta107,max_Omega106,max_Omega107\n"
                   <<std::scientific<<std::setprecision(12);
         csv<<"time_s,U_J,Psrc_W,Pperf_W,energy_balance_W,balance_rel_pct,Tmin_C,Tmax_C,CG_iter,CG_relres\n"
            <<std::scientific<<std::setprecision(12);
      }
      save(0,0.0);
      if(rank==0)
      {
         std::cout<<"Arrhenius kinetics: 106 A="<<in.A106<<" 1/s Ea="<<in.Ea106
                  <<" J/mol; 107 A="<<in.A107<<" 1/s Ea="<<in.Ea107<<" J/mol\n"
                  <<"Damage coverage threshold theta="<<in.damage_threshold
                  <<" (Omega >= "<<omega_cut<<")\n";
      }
      const int steps=static_cast<int>(std::llround(in.tend/in.dt));
      double Uprev=0;
      if(rank==0)
      {
         std::cout<<"Perfusion 106 [1/s]="<<in.omega106
                  <<"; 107 [1/s]="<<in.omega107
                  <<"; T0=Tb=37 C; adiabatic all thermal boundaries\n"
                  <<"time [s]   Tmax [C]   Pperf [W]  energy residual [%]\n";
      }
      for(int n=1;n<=steps;n++)
      {
         M->Mult(x,rhs); // P^T M P u_n (true DOF vector)
         rhs.Add(in.dt,*F); // RHS = M u_n + dt F_Q; T_b reference assumed 37C
         cg.Mult(rhs,x);
         const int converged_local=cg.GetConverged() ? 1 : 0;
         int converged=0;
         MPI_Allreduce(&converged_local,&converged,1,MPI_INT,MPI_MIN,comm);
         if(!converged)throw std::runtime_error("CG did not converge at heat step");
         rise.SetFromTrueDofs(x);
         temperature=rise;
         temperature+=in.initC;
         // These are exact FE integrals for spatially uniform heat capacity
         // and perfusion on each material, with linear nodal temperature.
         const double U=GlobalDot(comm,eform,rise);
         const double Pperf=GlobalDot(comm,pform,rise);
         const double imbalance=(U-Uprev)/in.dt + Pperf - Psrc;
         const double relative=100.0*imbalance/Psrc;
         Uprev=U;
         double maxloc=-std::numeric_limits<double>::infinity();
         double minloc= std::numeric_limits<double>::infinity();
         for(int i=0;i<temperature.Size();i++)
         {
            maxloc=std::max(maxloc,temperature(i));
            minloc=std::min(minloc,temperature(i));
         }
         const double tmin=GlobalMin(comm,minloc),tmax=GlobalMax(comm,maxloc);
         if(!std::isfinite(U)||!std::isfinite(tmax)||!std::isfinite(relative))
            throw std::runtime_error("Nonfinite thermal result");
         DamageStats d=AdvanceDamage(mesh,damageIR,rise_prev,rise,in.initC,in.dt,
                                    in.A106,in.Ea106,in.A107,in.Ea107,omega_cut,
                                    omega_at_qp,qfes,dgOmega,dgTheta,dgCoverage);
         ReduceDamage(comm,d);
         if(!std::isfinite(d.V107) || d.V107<0 || d.V107>sums.volume107*(1+1e-7) ||
            !std::isfinite(d.V106) || d.V106<0 || d.V106>sums.volume106*(1+1e-7))
            throw std::runtime_error("Invalid Arrhenius damaged volume");
         const bool checkpoint=(n%in.output_every==0 || n==steps);
         if(checkpoint)
         {
            const double time=n*in.dt;
            save(n,time);
            if(rank==0)
            {
               csv<<time<<','<<U<<','<<Psrc<<','<<Pperf<<','<<imbalance<<','
                  <<relative<<','<<tmin<<','<<tmax<<','<<cg.GetNumIterations()
                  <<','<<cg.GetFinalRelNorm()<<'\n';
               csv.flush();
               damage_csv<<time<<','<<in.damage_threshold<<','<<omega_cut<<','
                         <<1e6*d.V106<<','<<1e6*d.V107<<','
                         <<100.*d.V107/sums.volume107<<','
                         <<100.*d.V106/sums.volume106<<','<<1e6*d.V106<<','
                         <<d.thetaInt106/sums.volume106<<','
                         <<d.thetaInt107/sums.volume107<<','
                         <<d.maxOmega106<<','<<d.maxOmega107<<'\n';
               damage_csv.flush();
               std::cout<<std::fixed<<std::setprecision(3)
                        <<time<<"  "<<tmax<<"  "<<Pperf<<"  "<<relative
                        <<"  | tumor coverage: "<<100.*d.V107/sums.volume107
                        <<"% ; healthy damage: "<<1e6*d.V106<<" cm3\n";
            }
         }
         // A strict energy audit is made at every time step. A CG tolerance
         // affects this identity; report, not silently reinterpret an error.
         if(std::abs(relative)>0.1)
            throw std::runtime_error("Heat energy balance >0.1%; investigate FEM/solver");
         rise_prev=rise; // old temperature for next trapezoidal Arrhenius step
      }
      if(rank==0)
         std::cout<<"Done. Result: rezultati/"<<outName<<"/"<<outName<<".pvd\n"
                  <<"Energy: "<<outName<<"_energy.csv\n"
                  <<"Damage: "<<outName<<"_damage.csv\n";
   }
   catch(const std::exception &e)
   {
      std::cerr<<"[rank "<<rank<<"] ERROR: "<<e.what()<<'\n';
      MPI_Abort(comm,2);
      return 2;
   }
   MPI_Finalize();
   return 0;
}
