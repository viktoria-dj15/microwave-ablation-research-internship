// MPI verzija kontrolnog elektromagnetnog modela antene SA TUMOROM (107).
// Bazirano na korisnickom antena.cpp i MFEM v4.10 ex22p.cpp.
//
// Vazno: zadrzava iste fizicke koeficijente i granicne uslove kao
// serijska verzija. Izvozi Q za budući bioheat proračun, ali ne izracunava temperaturu ni ostecenje.
// Resavanje: paralelna direktna faktorizacija MUMPS.
// Ako provera originalnog sistema ne uspe, polje se NE snima.

#include "mfem.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <string>

#ifndef MFEM_USE_MPI
#error "antena_em_p.cpp zahteva MFEM izgradjen sa MFEM_USE_MPI=YES."
#endif
// U MFEM 4.10, HYPRE je zavisnost paralelne (MFEM_USE_MPI) izgradnje.
// Ne proveravati MFEM_USE_HYPRE: taj uslov nije pouzdan MFEM makro.
#ifndef MFEM_USE_MUMPS
#error "Ovaj program zahteva MFEM izgradjen sa MFEM_USE_MUMPS=ON."
#endif

namespace
{
constexpr int metal = 101;
constexpr int izolator = 102;
constexpr int slot1 = 103;
constexpr int slot2 = 104;
constexpr int kateter = 105;
constexpr int tkivo = 106;
constexpr int tumor = 107;
constexpr int port = 201;
constexpr int pec = 202;
constexpr int spolja = 203;
constexpr int kateter_ulaz = 204;
constexpr double coax_r_inner=0.25e-3;
constexpr double coax_r_outer=0.8e-3;
constexpr double coax_epsr=2.03;
constexpr double excitation_voltage=1.0; // Vpeak incidentnog TEM talasa; nije 1 W
constexpr double target_incident_power_W = 10.0; // nominalna INCIDENTNA snaga, nije neto snaga
constexpr mfem::ComplexOperator::Convention conv = mfem::ComplexOperator::HERMITIAN;

// TEM profil na ulazu koaksijalnog voda u ravni x=0, vršni napon 1 V.
// Ova normalizacija NE predstavlja ulaznu snagu 1 W.
void PortE(const mfem::Vector &x, mfem::Vector &e)
{
   constexpr double r_inner = coax_r_inner;
   constexpr double r_outer = coax_r_outer;
   constexpr double voltage = excitation_voltage;
   e.SetSize(3);
   e = 0.0;
   const double r2 = x[1] * x[1] + x[2] * x[2];
   if (r2 < r_inner * r_inner * 0.25) { return; }
   const double factor = voltage / (std::log(r_outer / r_inner) * r2);
   e[1] = factor * x[1];
   e[2] = factor * x[2];
}

// Imaginarna komponenta koaksijalne Robin-pobude, u obliku
// f = -n x (2k/mu E_inc), kako bi n x f = 2k/mu E_inc.
// Na portu 201 spoljna normala je n=(-1,0,0).
// VectorFEBoundaryTangentLFIntegrator integrise (n x f) . test.
double tem_forcing_k_over_mu=0.0;
void PortTEMImagForcing(const mfem::Vector &x, mfem::Vector &f)
{
   mfem::Vector e(3);
   PortE(x,e);
   f.SetSize(3);
   f=0.0;
   const double c=2.0*tem_forcing_k_over_mu;
   // -n x E = +ex x E za n = -ex.
   f[1]=-c*e[2];
   f[2]= c*e[1];
}

bool HasAttribute(const mfem::Array<int> &attributes, int attr)
{
   for (int i = 0; i < attributes.Size(); ++i)
   {
      if (attributes[i] == attr) { return true; }
   }
   return false;
}

int CountElements(const mfem::Mesh &mesh, int attr, bool boundary)
{
   const int count = boundary ? mesh.GetNBE() : mesh.GetNE();
   int selected = 0;
   for (int i = 0; i < count; ++i)
   {
      const int current = boundary ? mesh.GetBdrElement(i)->GetAttribute()
                                   : mesh.GetElement(i)->GetAttribute();
      selected += (current == attr);
   }
   return selected;
}

void DescribeBoundary(const mfem::Mesh &mesh, int attr)
{
   double xmin = std::numeric_limits<double>::infinity();
   double xmax = -xmin;
   double rmin = xmin;
   double rmax = 0.0;
   std::set<int> susedni_materijali;
   mfem::Array<int> vertices;
   for (int b = 0; b < mesh.GetNBE(); ++b)
   {
      if (mesh.GetBdrAttribute(b) != attr) { continue; }
      int element, info;
      mesh.GetBdrElementAdjacentElement(b, element, info);
      susedni_materijali.insert(mesh.GetAttribute(element));
      mesh.GetBdrElementVertices(b, vertices);
      for (int j = 0; j < vertices.Size(); ++j)
      {
         const mfem::real_t *p = mesh.GetVertex(vertices[j]);
         const double x = p[0];
         const double r = std::sqrt(p[1] * p[1] + p[2] * p[2]);
         xmin = std::min(xmin, x);
         xmax = std::max(xmax, x);
         rmin = std::min(rmin, r);
         rmax = std::max(rmax, r);
      }
   }
   std::cout << "  x temena [mm]: " << 1e3 * xmin << " do " << 1e3 * xmax
             << "; r temena [mm]: " << 1e3 * rmin << " do " << 1e3 * rmax
             << "\n  Susedne zapremine:";
   for (int material : susedni_materijali) { std::cout << ' ' << material; }
   std::cout << '\n';
}

// Proverava osnovnu geometriju otvorenih granica serialne EM podmreze.
bool CheckOpenBoundary(const mfem::Mesh &mesh, int attr, int material)
{
   mfem::Array<int> vertices;
   for (int b = 0; b < mesh.GetNBE(); ++b)
   {
      if (mesh.GetBdrAttribute(b) != attr) { continue; }
      int element, info;
      mesh.GetBdrElementAdjacentElement(b, element, info);
      if (mesh.GetAttribute(element) != material) { return false; }
      if (attr != kateter_ulaz) { continue; }
      mesh.GetBdrElementVertices(b, vertices);
      for (int j = 0; j < vertices.Size(); ++j)
      {
         const mfem::real_t *p = mesh.GetVertex(vertices[j]);
         const double r = std::sqrt(p[1] * p[1] + p[2] * p[2]);
         const double tolerance = 1e-7;
         if (std::abs(p[0]) > tolerance || r < 1e-3 - tolerance ||
             r > 1.2e-3 + tolerance) { return false; }
      }
   }
   return CountElements(mesh, attr, true) > 0;
}

// Validacija da je TEM port stvarno ulazni anulus u dielektriku 102.
// Ovo je obavezno za projekciju Z0, incidentnog i reflektovanog talasa.
bool CheckTEMPortGeometry(const mfem::Mesh &mesh)
{
   mfem::Array<int> v;
   int count=0;
   constexpr double tol=5e-7; // tolerancija CAD/geometrija u metrima
   for(int b=0;b<mesh.GetNBE();++b)
   {
      if(mesh.GetBdrAttribute(b)!=port)continue;
      ++count;
      int elem=-1,info=-1;
      mesh.GetBdrElementAdjacentElement(b,elem,info);
      if(elem<0 || mesh.GetAttribute(elem)!=izolator)return false;
      mesh.GetBdrElementVertices(b,v);
      for(int j=0;j<v.Size();++j)
      {
         const mfem::real_t *x=mesh.GetVertex(v[j]);
         const double r=std::sqrt(x[1]*x[1]+x[2]*x[2]);
         if(std::abs(x[0])>tol || r<coax_r_inner-tol ||
            r>coax_r_outer+tol)return false;
      }
   }
   return count>0;
}

// Norma distribuisanog vektora: nikad ne zamenjivati lokalnim Norml2().
double ParallelNorm(MPI_Comm comm, const mfem::Vector &x)
{
   const double local_norm = x.Norml2();
   const double local_square = local_norm * local_norm;
   double sum = 0.0;
   MPI_Allreduce(&local_square, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
   return std::sqrt(sum);
}

double ParallelMax(MPI_Comm comm, double local)
{
   double global = 0.0;
   MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, comm);
   return global;
}

bool ParallelAll(MPI_Comm comm, bool local)
{
   const int local_int = local ? 1 : 0;
   int global_int = 0;
   MPI_Allreduce(&local_int, &global_int, 1, MPI_INT, MPI_MIN, comm);
   return global_int != 0;
}

// Provere se obavljaju pre deljenja podmreze, da oznake budu identicne
// izvornoj serijskoj verziji. Sve MPI instance citaju isti .msh fajl.
bool ValidateSerialMesh(mfem::Mesh &cela_mreza, mfem::SubMesh &mesh, int rank)
{
   if (rank != 0) { return true; }
   const int uklonjeni_metal = CountElements(cela_mreza, metal, false);
   std::cout << "Cela mreza: " << cela_mreza.GetNE() << " elemenata\n"
             << "EM podmreza: " << mesh.GetNE() << " elemenata\n"
             << "Izostavljeni metal 101: " << uklonjeni_metal << " elemenata\n";
   if (mesh.GetNE() != cela_mreza.GetNE() - uklonjeni_metal)
   {
      std::cerr << "Neocekivan broj elemenata u EM podmrezi.\n";
      return false;
   }

   bool all_present = true;
   for (int attr : {izolator, slot1, slot2, kateter, tkivo, tumor})
   {
      const int count = CountElements(mesh, attr, false);
      std::cout << "Zapremina " << attr << ": " << count << '\n';
      all_present &= HasAttribute(mesh.attributes, attr) && count > 0;
   }
   if (CountElements(mesh, metal, false) != 0)
   {
      std::cerr << "Metal 101 ne sme biti u EM podmrezi.\n";
      return false;
   }
   for (int attr : {port, pec, spolja})
   {
      const int count = CountElements(mesh, attr, true);
      std::cout << "Povrsina " << attr << ": " << count << '\n';
      all_present &= HasAttribute(mesh.bdr_attributes, attr) && count > 0;
   }
   if (!all_present)
   {
      std::cerr << "Nedostaju ocekivani fizicki atributi.\n";
      return false;
   }
   for (int i = 0; i < mesh.bdr_attributes.Size(); ++i)
   {
      const int attr = mesh.bdr_attributes[i];
      if (attr != port && attr != pec && attr != spolja && attr != kateter_ulaz)
      {
         std::cerr << "Nepoznata granica " << attr << ".\n";
         DescribeBoundary(mesh, attr);
         return false;
      }
   }
   if(!CheckTEMPortGeometry(mesh))
   {
      std::cerr<<"Port 201 nije ispravan koaksijalni anulus na izolatoru 102.\n";
      return false;
   }
   if (!CheckOpenBoundary(mesh, spolja, tkivo) ||
       !CheckOpenBoundary(mesh, kateter_ulaz, kateter))
   {
      std::cerr << "Granice 203/204 nisu na ocekivanom materijalu.\n";
      return false;
   }
   std::cout << "Granica 204 (kateter): "
             << CountElements(mesh, kateter_ulaz, true) << " elemenata\n";
   DescribeBoundary(mesh, kateter_ulaz);
   return true;
}

int SaveField(mfem::ParMesh &mesh, mfem::ParComplexGridFunction &e,
              mfem::ParGridFunction &q_reference,
              mfem::ParGridFunction *q_incident_10W,
              mfem::ParGridFunction *er_incident_10W,
              mfem::ParGridFunction *ei_incident_10W,
              MPI_Comm comm, int rank, bool tem_port)
{
   // Za razliku od serijskog koda, nijedan proces ne sme da prepise
   // datoteke drugog procesa. ParaViewDataCollection to resava.
   const std::string output_name=tem_port?
      "antena_tumor_Pinc10W_Q_mpi":"antena_tumor_polje_mpi_fixed_port_Q1V";
   mfem::ParaViewDataCollection paraview(output_name, &mesh);
   paraview.SetPrefixPath("rezultati");
   paraview.SetLevelsOfDetail(1);
   paraview.SetDataFormat(mfem::VTKFormat::BINARY);
   paraview.SetCycle(0);
   paraview.SetTime(0.0); // stacionarni fazorski EM proracun
   paraview.RegisterField("E_real", &e.real());
   paraview.RegisterField("E_imag", &e.imag());
   // DG0 (L2 red 0): Q je konstantan po tetraedru i tacno cuva
   // zapreminski integral Q iz originalnog ND polja za ovaj red kvadrature.
   paraview.RegisterField("Q_ref_1Vpeak_W_m3", &q_reference);
   if (tem_port && q_incident_10W && er_incident_10W && ei_incident_10W)
   {
      paraview.RegisterField("Q_Pinc_10W_W_m3", q_incident_10W);
      paraview.RegisterField("E_real_Pinc_10W", er_incident_10W);
      paraview.RegisterField("E_imag_Pinc_10W", ei_incident_10W);
   }
   paraview.Save();
   const bool saved = ParallelAll(comm, paraview.Error() == 0);
   if (rank == 0)
   {
      if (saved)
      {
         std::cout << "Sacuvan ParaView izlaz: rezultati/" << output_name << "/.\n"
                   << (tem_port?"Incidentni TEM profil je zadat sa 1 Vpeak.\n":
                                  "Dirichlet ukupni E profil je zadat sa 1 Vpeak.\n");
      }
      else
      {
         std::cerr << "Greska pri paralelnom zapisivanju ParaView izlaza.\n";
      }
   }
   return saved ? 0 : 3;
}

// Rezultat zapreminske integracije po originalnim ND konacnim elementima.
// Nema projektovanja E u nodalno H1 polje, niti interpolacije preko ParaView.
struct PowerIntegrals
{
   double p106 = 0.0;
   double p107 = 0.0;
   double v106 = 0.0;
   double v107 = 0.0;
   double qmax106 = 0.0;
   double qmax107 = 0.0;
   long long n106 = 0;
   long long n107 = 0;
};

// Za svaki lokalni tetraedar: Q = 0.5 * sigma * (|Re E|^2 + |Im E|^2).
// Integrise Q * |det J| po referentnom tetraedru i sabira sve MPI procese.
// Integraciono pravilo reda 4 je iznad potrebnog reda za |E|^2 na
// afino mapiranom ND elementu prvog reda (vazi za nasu linearnu mrezu).
bool IntegrateTissuePower(mfem::ParMesh &mesh,
                          const mfem::ParComplexGridFunction &e,
                          double sigma106, double sigma107,
                          MPI_Comm comm, PowerIntegrals &result)
{
   std::array<double, 4> local_sum = {0.0, 0.0, 0.0, 0.0};
   std::array<double, 4> global_sum = {0.0, 0.0, 0.0, 0.0};
   std::array<double, 2> local_max = {0.0, 0.0};
   std::array<double, 2> global_max = {0.0, 0.0};
   long long local_count[2] = {0, 0};
   long long global_count[2] = {0, 0};
   bool local_ok = true;
   mfem::Vector er(mesh.SpaceDimension()), ei(mesh.SpaceDimension());

   for (int elem = 0; elem < mesh.GetNE(); ++elem)
   {
      const int attr = mesh.GetAttribute(elem);
      if (attr != tkivo && attr != tumor) { continue; }
      const int idx = attr == tkivo ? 0 : 1;
      const double sigma = idx == 0 ? sigma106 : sigma107;
      local_count[idx]++;

      mfem::ElementTransformation *trans = mesh.GetElementTransformation(elem);
      const mfem::IntegrationRule &ir =
         mfem::IntRules.Get(mesh.GetElementGeometry(elem), 4);
      for (int j = 0; j < ir.GetNPoints(); ++j)
      {
         const mfem::IntegrationPoint &ip = ir.IntPoint(j);
         trans->SetIntPoint(&ip);

         // Evalucija H(curl) polja U UNUTRASNJOSTI elementa (ND/Piola).
         // Ovaj overload je bezbedan u integracionoj petlji sa istim trans.
         e.real().GetVectorValue(*trans, ip, er);
         e.imag().GetVectorValue(*trans, ip, ei);
         if (er.Size() != mesh.SpaceDimension() ||
             ei.Size() != mesh.SpaceDimension())
         {
            local_ok = false;
            break;
         }
         double e2 = 0.0;
         for (int k = 0; k < er.Size(); ++k)
         {
            e2 += er[k] * er[k] + ei[k] * ei[k];
         }
         const double dv = ip.weight * trans->Weight(); // m^3
         const double q  = 0.5 * sigma * e2;             // W/m^3
         if (!std::isfinite(q) || !std::isfinite(dv) || dv <= 0.0)
         {
            local_ok = false;
            break;
         }
         local_sum[idx]     += q * dv; // W
         local_sum[idx + 2] += dv;     // m^3
         local_max[idx] = std::max(local_max[idx], q);
      }
      if (!local_ok) { break; }
   }

   // SVE instance moraju ucestvovati u MPI operacijama, cak i pri gresci.
   if (!ParallelAll(comm, local_ok)) { return false; }
   MPI_Allreduce(local_sum.data(), global_sum.data(), 4, MPI_DOUBLE,
                 MPI_SUM, comm);
   MPI_Allreduce(local_max.data(), global_max.data(), 2, MPI_DOUBLE,
                 MPI_MAX, comm);
   MPI_Allreduce(local_count, global_count, 2, MPI_LONG_LONG,
                 MPI_SUM, comm);

   result.p106 = global_sum[0];
   result.p107 = global_sum[1];
   result.v106 = global_sum[2];
   result.v107 = global_sum[3];
   result.qmax106 = global_max[0];
   result.qmax107 = global_max[1];
   result.n106 = global_count[0];
   result.n107 = global_count[1];
   return result.n106 > 0 && result.n107 > 0 &&
          result.v106 > 0.0 && result.v107 > 0.0 &&
          std::isfinite(result.p106) && std::isfinite(result.p107);
}

// ================= NOVO: IZVOR TOPLOTE ZA BIOHEAT =================
// Ne resava bioheat jednacinu. Q ostaje u W/m^3 i obuhvata samo materijale
// 106 i 107. Ostali EM materijali (izolator, slotovi, kateter) imaju Q=0,
// jer je njihova elektricna provodnost u ovoj simulaciji nula.
//
// Konstrukcija DG0/L2(0) polja: za svaki tetraedar se racuna zapreminski
// prosek q = integral(Q dV)/integral(dV), ne vrednost u centroidu.
// Time je integral iz ParaView DG0 izvoza konzervativan po elementima.
// Q ipak nije tackasta originalna funkcija Q(x): za termalni FEM koji koristi
// istu mrezu idealno je direktno evaluirati ND E u kvadraturnim tackama.
struct ThermalSourceVerification
{
   double p106_ref = 0.0;
   double p107_ref = 0.0;
   double p106_10W = 0.0;
   double p107_10W = 0.0;
   double qmax_cell_10W = 0.0;
};

bool MakeThermalSourceDG0(mfem::ParMesh &mesh,
                          const mfem::ParComplexGridFunction &e,
                          mfem::ParFiniteElementSpace &q_fes,
                          mfem::ParGridFunction &q_ref,
                          mfem::ParGridFunction *q_10W,
                          double sigma106, double sigma107,
                          double q_scale,
                          MPI_Comm comm,
                          ThermalSourceVerification &report)
{
   q_ref = 0.0;
   if (q_10W) { *q_10W = 0.0; }
   std::array<double,4> local = {{0.0,0.0,0.0,0.0}};
   std::array<double,4> global = {{0.0,0.0,0.0,0.0}};
   double qmax_local = 0.0, qmax_global = 0.0;
   bool ok = true;
   mfem::Array<int> edofs;
   mfem::Vector er(mesh.SpaceDimension()), ei(mesh.SpaceDimension());
   for (int elem=0; elem<mesh.GetNE(); ++elem)
   {
      q_fes.GetElementDofs(elem, edofs);
      if (edofs.Size()!=1 || edofs[0]<0)
      {
         ok=false;
         break;
      }
      const int attribute=mesh.GetAttribute(elem);
      const bool is_tissue=(attribute==tkivo || attribute==tumor);
      if (!is_tissue) { continue; } // ostaje nula (sigma=0)
      const double sigma=attribute==tkivo ? sigma106 : sigma107;
      auto *tr=mesh.GetElementTransformation(elem);
      const mfem::IntegrationRule &ir=mfem::IntRules.Get(mesh.GetElementGeometry(elem),4);
      double volume=0.0, power=0.0;
      for (int j=0; j<ir.GetNPoints(); ++j)
      {
         const auto &ip=ir.IntPoint(j);
         tr->SetIntPoint(&ip);
         e.real().GetVectorValue(*tr,ip,er);
         e.imag().GetVectorValue(*tr,ip,ei);
         double norm2=0.0;
         for (int k=0;k<er.Size();++k)
            norm2+=er[k]*er[k]+ei[k]*ei[k];
         const double dv=ip.weight*tr->Weight();
         const double q=0.5*sigma*norm2;
         if (!(dv>0.0) || !std::isfinite(dv) ||
             !std::isfinite(q) || q<0.0)
         {
            ok=false;
            break;
         }
         power+=q*dv;
         volume+=dv;
      }
      if (!ok || !(volume>0.0)) {ok=false;break;}
      const double q_mean=power/volume;
      q_ref(edofs[0])=q_mean;
      if (q_10W)
      {
         (*q_10W)(edofs[0])=q_mean*q_scale;
         qmax_local=std::max(qmax_local,q_mean*q_scale);
      }
      const int idx=(attribute==tkivo) ? 0 : 1;
      local[idx]+=q_mean*volume;
      if (q_10W) { local[idx+2]+=q_mean*q_scale*volume; }
   }
   if (!ParallelAll(comm,ok)) { return false; }
   MPI_Allreduce(local.data(),global.data(),4,MPI_DOUBLE,MPI_SUM,comm);
   MPI_Allreduce(&qmax_local,&qmax_global,1,MPI_DOUBLE,MPI_MAX,comm);
   report.p106_ref=global[0];
   report.p107_ref=global[1];
   report.p106_10W=global[2];
   report.p107_10W=global[3];
   report.qmax_cell_10W=qmax_global;
   return std::isfinite(report.p106_ref) && std::isfinite(report.p107_ref) &&
          std::isfinite(report.p106_10W) && std::isfinite(report.p107_10W);
}

bool WriteNominalPower10WCSV(const PowerIntegrals &p,
                             const ThermalSourceVerification &s,
                             const double z0, const double p_inc_nominal,
                             const double e_scale, const double q_scale,
                             const bool tem_port)
{
   if (!tem_port) {return true;}
   std::ofstream csv("fem_izvor_toplote_Pinc10W.csv");
   if (!csv) {return false;}
   csv << std::scientific << std::setprecision(14);
   csv << "kljuc,vrednost,jedinica\n"
       << "tip_normalizacije,nominalna_incidentna_TEM,-\n"
       << "zadati_Vinc_peak,"<<excitation_voltage<<",V\n"
       << "Z0,"<<z0<<",ohm\n"
       << "Pinc_nominal_ref,"<<p_inc_nominal<<",W\n"
       << "Pinc_cilj,"<<target_incident_power_W<<",W\n"
       << "faktor_E,"<<e_scale<<",1\n"
       << "faktor_Q,"<<q_scale<<",1\n"
       << "P106_ref_ND,"<<p.p106<<",W\n"
       << "P107_ref_ND,"<<p.p107<<",W\n"
       << "P106_DG0_ref,"<<s.p106_ref<<",W\n"
       << "P107_DG0_ref,"<<s.p107_ref<<",W\n"
       << "P106_10W,"<<s.p106_10W<<",W\n"
       << "P107_10W,"<<s.p107_10W<<",W\n"
       << "P_aps_10W,"<<(s.p106_10W+s.p107_10W)<<",W\n"
       << "Qmax_DG0_10W,"<<s.qmax_cell_10W<<",W_per_m3\n"
       << "Qmax_kvadratura106_10W,"<<(p.qmax106*q_scale)<<",W_per_m3\n"
       << "Qmax_kvadratura107_10W,"<<(p.qmax107*q_scale)<<",W_per_m3\n";
   return static_cast<bool>(csv);
}

// Pise CSV samo rank 0; korisno za izvestaj i poredjenje sa ParaView.
bool SavePowerCSV(const PowerIntegrals &p,
                  double sigma106, double sigma107, bool tem_port)
{
   std::ofstream out(tem_port?
      "fem_apsorbovana_snaga_106_107_TEM.csv":
      "fem_apsorbovana_snaga_106_107_fixed.csv");
   if (!out) { return false; }
   out << std::scientific << std::setprecision(12);
   out << "materijal,atribut,sigma_S_m,zapremina_m3,broj_tetraedara,"
          "P_aps_W,Q_prosek_W_m3,Q_maks_kvadratura_W_m3\n";
   out << "Zdravo_tkivo,106," << sigma106 << ',' << p.v106 << ','
       << p.n106 << ',' << p.p106 << ',' << p.p106 / p.v106 << ','
       << p.qmax106 << '\n';
   out << "Tumor,107," << sigma107 << ',' << p.v107 << ','
       << p.n107 << ',' << p.p107 << ',' << p.p107 / p.v107 << ','
       << p.qmax107 << '\n';
   return static_cast<bool>(out);
}

void PrintPower(const PowerIntegrals &p, int rank, bool tem_port)
{
   if (rank != 0) { return; }
   const double sum = p.p106 + p.p107;
   std::cout << std::scientific << std::setprecision(10)
             << "\n=== DIREKTNA MFEM INTEGRACIJA APSORBOVANE SNAGE ===\n"
             << "  Q = 0.5 * sigma * (|Re E|^2 + |Im E|^2)\n"
             << "  Normalizacija: "
             << (tem_port?"incidentni TEM 1 Vpeak":"ukupni E Dirichlet 1 Vpeak")
             << "; NISU terapijski vati.\n"
             << "  [106] Zdravo tkivo:\n"
             << "    Broj tetraedara: " << p.n106 << '\n'
             << "    Zapremina [mm^3]: " << p.v106*1e9 << '\n'
             << "    P_abs [W]: " << p.p106 << '\n'
             << "    Q_prosek [W/m^3]: " << p.p106/p.v106 << '\n'
             << "    Q_maks na kvadraturnim tackama [W/m^3]: " << p.qmax106 << '\n'
             << "  [107] Tumor:\n"
             << "    Broj tetraedara: " << p.n107 << '\n'
             << "    Zapremina [mm^3]: " << p.v107*1e9 << '\n'
             << "    P_abs [W]: " << p.p107 << '\n'
             << "    Q_prosek [W/m^3]: " << p.p107/p.v107 << '\n'
             << "    Q_maks na kvadraturnim tackama [W/m^3]: " << p.qmax107 << '\n'
             << "  P_abs_106+107 [W]: " << sum << '\n'
             << "  Udeo P_abs u tumoru [%]: " << (sum > 0 ? 100*p.p107/sum : 0.0)
             << '\n'
             << "  NAPOMENA: Zbir ne predstavlja automatski ulaznu snagu porta.\n"
             << "===================================================\n";
}


// ---------------------------------------------------------------------------
// Granicni Poyntingovi tokovi preko PRAVOG H(curl) / Nedelec polja.
// Fizicka fazorska konvencija operatora je exp(+i*omega*t):
//     curl E = -i*omega*mu*H,  H = i/(omega*mu) curl E.
//     H_re = -curl(E_im)/(omega*mu),
//     H_im = +curl(E_re)/(omega*mu).
// <S> = 0.5*Re(E x conj(H)) = 0.5*(E_re x H_re + E_im x H_im).
// Na granicnim trouglovima normalu orijentisemo SPOLJA iz volumena
// prema polozaju centra susednog tetraedra. Vektor normale nosi dA.
// Integraciono pravilo reda 4 koristi se za ND red 1 na affine mrezi.
struct BoundaryPowerFlux
{
   // Indeksi 0=port201, 1=PEC202, 2=ABC203, 3=ABC204.
   std::array<double,4> out = {{0.,0.,0.,0.}};
   std::array<double,4> area = {{0.,0.,0.,0.}};
   std::array<long long,4> faces = {{0,0,0,0}};
   // Na ABC, nezavisna evaluacija izlazne snage iz tangencijalnog E
   // i lokalne TE/plane-wave impedanse (prvi red ABC).
   double abc_e203=0.0,abc_e204=0.0;
   // TEM modalna projekcija na geometrijski diskretizovan port 201.
   std::complex<double> ve={0.,0.}, vh={0.,0.};
   std::complex<double> vinc={0.,0.}, vref={0.,0.};
   std::complex<double> s11={0.,0.};
   double mode_norm=0., modal_inc=0.,modal_ref=0.,modal_net=0.;
   double z0=0.;
   // Port L2 norme na povrsini 201; komponenta normalna na presek se izostavlja.
   double port_et_l2=0.,port_eta_ht_l2=0.;
};

int BoundaryFluxIndex(int attr)
{
   switch (attr)
   {
      case port:return 0;
      case pec:return 1;
      case spolja:return 2;
      case kateter_ulaz:return 3;
      default:return -1;
   }
}

// Metod A: evaluacija E i curl(E) na volumenskoj strani svakog granicnog
// tetraedra koristeci preslikavanje referentni-trougao -> tetraedar.
// Metod B: prethodna implementacija koja koristi boundary transformation.
// Poredimo A/B pri redovima kvadrature 2,4,6,8 (ista FEM solucija).
bool IntegratePoyntingFlux(mfem::ParMesh &mesh,
                           const mfem::ParComplexGridFunction &e,
                           double omega, double mu0, double epsilon0,
                           MPI_Comm comm, int quad_order, bool volume_side,
                           BoundaryPowerFlux &result)
{
   // 0..3 power, 4..7 area, 8..9 power via ABC(E),
   // 10..11 modal E, 12..13 modal eta*H, 14 modal norm,
   // 15 integral |E_t|^2, 16 integral |eta*H_t|^2.
   std::array<double,17> local = {{0.}};
   std::array<double,17> global = {{0.}};
   long long local_n[4]={0,0,0,0}, global_n[4]={0,0,0,0};
   bool local_ok=true;
   const double eta_coax=std::sqrt(mu0/(epsilon0*coax_epsr));
   const double logab=std::log(coax_r_outer/coax_r_inner);
   const std::complex<double> eps_203(epsilon0*57.55,-1.95/omega);
   const std::complex<double> eps_204(epsilon0*2.6,0.);
   const double y203=std::real(std::sqrt(eps_203/mu0));
   const double y204=std::real(std::sqrt(eps_204/mu0));
   mfem::Vector er(3),ei(3),cr(3),ci(3),hr(3),hi(3),normal(3),xyz(3);
   mfem::Array<int> verts;

   for (int b=0;b<mesh.GetNBE();++b)
   {
      int idx=BoundaryFluxIndex(mesh.GetBdrAttribute(b));
      if(idx<0) {local_ok=false;break;}
      int elem=-1,info=-1;
      mesh.GetBdrElementAdjacentElement(b,elem,info);
      if(elem<0 || elem>=mesh.GetNE()) {local_ok=false;break;}
      mesh.GetElementVertices(elem,verts);
      if(verts.Size()!=4) {local_ok=false;break;}
      double xc[3]={0.,0.,0.};
      for(int v=0;v<4;++v)
      {
         const mfem::real_t *xv=mesh.GetVertex(verts[v]);
         for(int k=0;k<3;++k) xc[k]+=0.25*xv[k];
      }

      mfem::FaceElementTransformations *ftr=nullptr;
      mfem::IsoparametricTransformation btr;
      mfem::ElementTransformation *face=nullptr;
      if(volume_side)
      {
         ftr=mesh.GetBdrFaceTransformations(b);
         if(!ftr || !ftr->Elem1 || !ftr->Face)
         {local_ok=false;break;}
         face=ftr->Face;
      }
      else
      {
         mesh.GetBdrElementTransformation(b,&btr);
         face=&btr;
      }
      const mfem::IntegrationRule &ir=
          mfem::IntRules.Get(mesh.GetBdrElementGeometry(b),quad_order);
      local_n[idx]++;
      for(int j=0;j<ir.GetNPoints();++j)
      {
         const mfem::IntegrationPoint &ip=ir.IntPoint(j);
         face->SetIntPoint(&ip);
         if(volume_side)
         {
            mfem::IntegrationPoint ip_elem;
            ftr->Loc1.Transform(ip,ip_elem);
            ftr->Elem1->SetIntPoint(&ip_elem);
            e.real().GetVectorValue(*ftr->Elem1,ip_elem,er);
            e.imag().GetVectorValue(*ftr->Elem1,ip_elem,ei);
            e.real().GetCurl(*ftr->Elem1,cr);
            e.imag().GetCurl(*ftr->Elem1,ci);
         }
         else
         {
            e.real().GetVectorValue(btr,ip,er);
            e.imag().GetVectorValue(btr,ip,ei);
            e.real().GetCurl(btr,cr);
            e.imag().GetCurl(btr,ci);
         }
         if(er.Size()!=3 || ei.Size()!=3 || cr.Size()!=3 || ci.Size()!=3)
         {local_ok=false;break;}
         for(int k=0;k<3;++k)
         {
            hr[k]=-ci[k]/(omega*mu0);
            hi[k]= cr[k]/(omega*mu0);
         }
         mfem::CalcOrtho(face->Jacobian(),normal); // normal * dA
         face->Transform(ip,xyz);
         double orientation=0.0;
         for(int k=0;k<3;++k) orientation+=normal[k]*(xyz[k]-xc[k]);
         if(!std::isfinite(orientation) || std::abs(orientation)<1e-24)
         {local_ok=false;break;}
         if(orientation<0.) normal*=-1.;

         const double sx=0.5*(er[1]*hr[2]-er[2]*hr[1]+ei[1]*hi[2]-ei[2]*hi[1]);
         const double sy=0.5*(er[2]*hr[0]-er[0]*hr[2]+ei[2]*hi[0]-ei[0]*hi[2]);
         const double sz=0.5*(er[0]*hr[1]-er[1]*hr[0]+ei[0]*hi[1]-ei[1]*hi[0]);
         const double dA=normal.Norml2()*ip.weight;
         const double dP=(sx*normal[0]+sy*normal[1]+sz*normal[2])*ip.weight;
         if(!std::isfinite(dP)||!std::isfinite(dA)||dA<=0.)
         {local_ok=false;break;}
         local[idx]+=dP;
         local[4+idx]+=dA;

         // Po prirodi prvorednog apsorbujuceg Robin uslova, granicna
         // disipacija se moze proceniti i samo iz tangencijalnog E.
         if(idx==2 || idx==3)
         {
            const double nn=normal.Norml2();
            double nt_r=0.,nt_i=0.,e2=0.;
            for(int k=0;k<3;++k)
            {
               nt_r+=er[k]*normal[k]/nn;
               nt_i+=ei[k]*normal[k]/nn;
               e2+=er[k]*er[k]+ei[k]*ei[k];
            }
            const double et2=std::max(0.,e2-nt_r*nt_r-nt_i*nt_i);
            local[idx==2?8:9]+=0.5*(idx==2?y203:y204)*et2*dA;
         }

         if(idx==0)
         {
            const double r2=xyz[1]*xyz[1]+xyz[2]*xyz[2];
            if(!(r2>0.0)) {local_ok=false;break;}
            // E_mode=er/(r ln(b/a)), norm integral 2pi/ln(b/a).
            const double ty=xyz[1]/(r2*logab);
            const double tz=xyz[2]/(r2*logab);
            const double phiy=-tz,phiz=ty;
            local[10]+=(er[1]*ty+er[2]*tz)*dA;
            local[11]+=(ei[1]*ty+ei[2]*tz)*dA;
            local[12]+=eta_coax*(hr[1]*phiy+hr[2]*phiz)*dA;
            local[13]+=eta_coax*(hi[1]*phiy+hi[2]*phiz)*dA;
            local[14]+=(ty*ty+tz*tz)*dA;
            // Lokalni poprecni intenziteti; korisno za kvantifikaciju
            // dela polja koji NE pripada osnovnom koaksijalnom TEM modu.
            double ntEre=0.,ntEim=0.,ntHre=0.,ntHim=0.;
            double e2=0.,h2=0.;
            const double invnormal=1./normal.Norml2();
            for(int k=0;k<3;++k)
            {
               const double n=normal[k]*invnormal;
               ntEre+=er[k]*n; ntEim+=ei[k]*n;
               ntHre+=hr[k]*n; ntHim+=hi[k]*n;
               e2+=er[k]*er[k]+ei[k]*ei[k];
               h2+=hr[k]*hr[k]+hi[k]*hi[k];
            }
            local[15]+=(e2-ntEre*ntEre-ntEim*ntEim)*dA;
            local[16]+=(h2-ntHre*ntHre-ntHim*ntHim)*eta_coax*eta_coax*dA;
         }
      }
      if(!local_ok) break;
   }
   if(!ParallelAll(comm,local_ok)) return false;
   MPI_Allreduce(local.data(),global.data(),17,MPI_DOUBLE,MPI_SUM,comm);
   MPI_Allreduce(local_n,global_n,4,MPI_LONG_LONG,MPI_SUM,comm);
   for(int k=0;k<4;++k)
   {
      result.out[k]=global[k];result.area[k]=global[4+k];result.faces[k]=global_n[k];
      if(result.faces[k]<=0||!std::isfinite(result.out[k])||result.area[k]<=0.)return false;
   }
   result.abc_e203=global[8];result.abc_e204=global[9];
   result.mode_norm=global[14];
   result.port_et_l2=global[15];
   result.port_eta_ht_l2=global[16];
   if(!(result.mode_norm>0.)||!std::isfinite(result.mode_norm))return false;
   result.ve={global[10]/result.mode_norm,global[11]/result.mode_norm};
   result.vh={global[12]/result.mode_norm,global[13]/result.mode_norm};
   result.vinc=0.5*(result.ve+result.vh);
   result.vref=0.5*(result.ve-result.vh);
   result.z0=eta_coax*logab/(2.*std::acos(-1.));
   result.modal_inc=std::norm(result.vinc)/(2.*result.z0);
   result.modal_ref=std::norm(result.vref)/(2.*result.z0);
   result.modal_net=result.modal_inc-result.modal_ref;
   if(std::abs(result.vinc)>1e-15)result.s11=result.vref/result.vinc;
   return std::isfinite(result.modal_inc)&&std::isfinite(result.modal_ref);
}


struct PowerBalance
{
   double pin_net = 0.0;
   double pabs = 0.0;
   double p203 = 0.0;
   double p204 = 0.0;
   double ppec = 0.0;
   double residual = 0.0;
   double relative = 0.0;
};

PowerBalance CalculatePowerBalance(const PowerIntegrals &p,
                                   const BoundaryPowerFlux &f)
{
   PowerBalance b;
   b.pin_net = -f.out[0];
   b.pabs = p.p106+p.p107;
   b.p203 = f.out[2];
   b.p204 = f.out[3];
   b.ppec = f.out[1];
   b.residual = b.pin_net - b.pabs - b.p203 - b.p204 - b.ppec;
   // Delilac nije samo P_in: ako je blizu nule, izbegavamo nestabilan
   // procenat. Ova mera opisuje numericku gresku bilansa, a ne efikasnost.
   const double scale = std::max({std::abs(b.pin_net), b.pabs,
         std::abs(b.p203)+std::abs(b.p204)+std::abs(b.ppec), 1e-30});
   b.relative = std::abs(b.residual)/scale;
   return b;
}

void PrintPowerBalance(const PowerIntegrals &p,
                       const BoundaryPowerFlux &f,
                       const PowerBalance &b, int rank, bool tem_port)
{
   if (rank != 0) { return; }
   constexpr double previous_absorption_W = 7.9478612835e-3;
   std::cout << std::scientific << std::setprecision(10)
      << "\n=== POYNTING / ELEKTROMAGNETNI ENERGETSKI BILANS ===\n"
      << "  Konvencija e^(+i*omega*t), H=i curl(E)/(omega*mu).\n"
      << "  Pozitivan fluks = izlaz IZ EM domena (spoljasnja normala).\n"
      << "  201 port : tok spolja pozitivno [W]: " << f.out[0] << '\n'
      << "  201 port : NETO ulaz = -tok201 [W]: " << b.pin_net << '\n'
      << "  203 ABC  : izlazni tok [W]: " << b.p203 << '\n'
      << "  204 ABC  : izlazni tok [W]: " << b.p204 << '\n'
      << "  202 PEC  : tok (treba ~0) [W]: " << b.ppec << '\n'
      << "  P_abs 106 [W]: " << p.p106 << '\n'
      << "  P_abs 107 [W]: " << p.p107 << '\n'
      << "  P_abs ukupno 106+107 [W]: " << b.pabs << '\n'
      << "  BILANS: P_neto201 - P_abs - P_iz203 - P_iz204 - P_202 [W]: "
      << b.residual << '\n'
      << "  Relativno apsolutno odstupanje bilansa [%]: "
      << b.relative*100.0 << '\n'
      << "  [Broj granicnih trouglova] 201=" << f.faces[0]
      << "  202=" << f.faces[1]
      << "  203=" << f.faces[2]
      << "  204=" << f.faces[3] << '\n'
      << "  Ovo je NETO snaga porta, NE incidentna snaga generatora.\n"
      << "  Bilans i mala algebrajska greska nisu dokaz konvergencije mreze/ABC.\n"
      << "==========================================================\n";
   if (!tem_port)
   {
      std::cout << "  Razlika od prethodne apsorpcije [W]: "
                << b.pabs-previous_absorption_W << '\n';
   }
   if (b.pin_net <= 0.0)
   {
      std::cout << "  UPOZORENJE: neto ulaz kroz port nije pozitivan;"
                   " proveriti znakove, granice i pobudu.\n";
   }
   if (b.p203 < -1e-9 || b.p204 < -1e-9)
   {
      std::cout << "  UPOZORENJE: negativan izlaz kroz ABC;"
                   " proveriti granicni uslov i refleksije.\n";
   }
   if (b.relative > 0.05)
   {
      std::cout << "  UPOZORENJE: razlika bilansa prelazi 5%."
                   " Pre termickog modela proveriti tokove, ABC i diskretizaciju.\n";
   }
}

bool SavePowerBalanceCSV(const PowerIntegrals &p,
                         const BoundaryPowerFlux &f,
                         const PowerBalance &b, bool tem_port)
{
   std::ofstream out(tem_port?
      "fem_poynting_bilans_TEM.csv":"fem_poynting_bilans_fixed.csv");
   if (!out) { return false; }
   out << std::scientific << std::setprecision(12);
   out << "velicina,vrednost_W,napomena\n";
   out << "fluks201_spolja," << f.out[0] << ",spoljasnja_normala\n";
   out << "port201_neto_ulaz," << b.pin_net << ",minus_fluks201\n";
   out << "fluks202_PEC," << b.ppec << ",idealno_nula\n";
   out << "fluks203_ABC," << b.p203 << ",izlaz_spolja\n";
   out << "fluks204_ABC," << b.p204 << ",izlaz_spolja\n";
   out << "apsorpcija106," << p.p106 << ",direktna_MFEM_kvadratura\n";
   out << "apsorpcija107," << p.p107 << ",direktna_MFEM_kvadratura\n";
   out << "apsorpcija106_107," << b.pabs << ",zbir\n";
   out << "rezidual_bilansa," << b.residual << ",201_apsorpcija_203_204_202\n";
   out << "relativni_rezidual_bilansa," << b.relative
       << ",bezdimenziono\n";
   return static_cast<bool>(out);
}


// ---------------------------------------------------------------------------
// Dijagnostika 1: diskretna slaba energija iz SISTEMA koji je MUMPS resio.
// Za TEM Robin (homogeni PEC): E pripada prostoru test funkcija. Vazi
//   Im( X^H A X ) / (2 omega)  =  Im( X^H B ) / (2 omega).
// Na levoj strani doprinose daju omega*sigma i Im(q_ABC), Im(q_port).
// Zato: P_weak_port = Im(X^H B)/(2omega) - P_port_Robin.
//
// Fixed port ima nenulti zadati Dirichlet podatak i E nije dopustena
// test funkcija. Za fixed se izracunava samo algebrajski identitet,
// bez tumacenja X^H B kao fizicke snage.
struct WeakEnergyDiagnostics
{
   std::complex<double> xhax={0.,0.},xhb={0.,0.};
   double power_im_ax=0.,power_im_b=0.,power_im_diff=0.;
   double p_robin201=0.;
   double p_abc203=0.,p_abc204=0.;
   double p_materials=0.;
   double p_weak_port=0.;
   double p_weak_vs_poynting=0.;
   double p_matrix_vs_integrals=0.;
};

bool CalculateWeakEnergy(const mfem::Operator &A,
                         const mfem::Vector &X,
                         const mfem::Vector &B,
                         double omega, double epsilon0, double mu0,
                         const BoundaryPowerFlux &f,
                         const PowerIntegrals &p,
                         MPI_Comm comm, bool tem_port,
                         WeakEnergyDiagnostics &d)
{
   mfem::Vector AX(X.Size());
   A.Mult(X, AX);
   if(X.Size()%2!=0 || B.Size()!=X.Size()) return false;
   const int n=X.Size()/2;
   double local[4]={0.,0.,0.,0.},glob[4]={0.,0.,0.,0.};
   for(int i=0;i<n;++i)
   {
      const double xr=X[i], xi=X[n+i];
      local[0]+=xr*AX[i]+xi*AX[n+i];
      local[1]+=xr*AX[n+i]-xi*AX[i];
      local[2]+=xr*B[i]+xi*B[n+i];
      local[3]+=xr*B[n+i]-xi*B[i];
   }
   MPI_Allreduce(local,glob,4,MPI_DOUBLE,MPI_SUM,comm);
   d.xhax={glob[0],glob[1]}; d.xhb={glob[2],glob[3]};
   d.power_im_ax=glob[1]/(2.*omega);
   d.power_im_b=glob[3]/(2.*omega);
   d.power_im_diff=d.power_im_ax-d.power_im_b;
   const double eta=std::sqrt(mu0/(epsilon0*coax_epsr));
   d.p_robin201=tem_port?f.port_et_l2/(2.*eta):0.;
   d.p_materials=p.p106+p.p107;
   d.p_abc203=f.abc_e203;
   d.p_abc204=f.abc_e204;
   // Ovaj identitet vazi samo za TEM, uz nulti tangencijalni E na PEC.
   if(tem_port)
   {
      d.p_weak_port=d.power_im_b-d.p_robin201;
      d.p_weak_vs_poynting=d.p_weak_port+f.out[0];
      d.p_matrix_vs_integrals=d.power_im_ax-
          (d.p_materials+d.p_abc203+d.p_abc204+d.p_robin201);
   }
   return std::isfinite(d.power_im_ax) && std::isfinite(d.power_im_b);
}

void WriteWeakEnergy(const WeakEnergyDiagnostics &d,
                     const BoundaryPowerFlux &f,
                     int rank, bool tem_port)
{
   if(rank!=0)return;
   const std::string suffix=tem_port?"TEM":"fixed";
   std::cout<<"\n=== DISKRETNA SLABA ENERGETSKA PROVERA (A,X,B) ===\n"
            <<"  Im(X^H A X)/(2omega) [W]: "<<d.power_im_ax<<'\n'
            <<"  Im(X^H B)/(2omega) [W]: "<<d.power_im_b<<'\n'
            <<"  Razlika leva/desna slaba forma [W]: "<<d.power_im_diff<<'\n';
   if(tem_port)
   {
      std::cout<<"  P_ABC203 iz slabog Robin izraza [W]: "<<d.p_abc203<<'\n'
               <<"  P_ABC204 iz slabog Robin izraza [W]: "<<d.p_abc204<<'\n'
               <<"  P_robin201 = int|Et|^2/(2 eta) [W]: "<<d.p_robin201<<'\n'
               <<"  P_tkiva+P_ABC+P_robin201 [W]: "
               <<d.p_materials+d.p_abc203+d.p_abc204+d.p_robin201<<'\n'
               <<"  Matrica naspram FEM integrala [W]: "
               <<d.p_matrix_vs_integrals<<'\n'
               <<"  P_neto201 po slaboj formulaciji [W]: "<<d.p_weak_port<<'\n'
               <<"  P_neto201 direktni Poynting [W]: "<<-f.out[0]<<'\n'
               <<"  Razlika weak-port minus Poynting [W]: "
               <<d.p_weak_vs_poynting<<'\n';
   }
   else
   {
      std::cout<<"  NENULTI Dirichlet port: Im(X^H B)/(2omega) NIJE "
                  "fizicka ulazna snaga. Samo algebrajski test.\n";
   }
   std::cout<<"====================================================\n";
   std::ofstream out("fem_slaba_energija_"+suffix+".csv");
   if(!out)return;
   out<<std::scientific<<std::setprecision(13)
      <<"rezim,Im_XHAX_po_2omega_W,Im_XHB_po_2omega_W,algebra_razlika_W,"
        "P_abs_W,P_abc203_E_W,P_abc204_E_W,P_robin_port201_W,"
        "P_weak_neto201_W,P_poynting_neto201_W,weak_poynting_razlika_W,"
        "A_vs_integrali_W,validan_weak_port\n"
      <<(tem_port?"TEM":"fixed")<<','<<d.power_im_ax<<','<<d.power_im_b<<','
      <<d.power_im_diff<<','<<d.p_materials<<','<<d.p_abc203<<','
      <<d.p_abc204<<','<<d.p_robin201<<','
      <<(tem_port?d.p_weak_port:std::numeric_limits<double>::quiet_NaN())<<','
      <<-f.out[0]<<','
      <<(tem_port?d.p_weak_vs_poynting:std::numeric_limits<double>::quiet_NaN())<<','
      <<(tem_port?d.p_matrix_vs_integrals:std::numeric_limits<double>::quiet_NaN())<<','
      <<(tem_port?1:0)<<'\n';
}

// ---------------------------------------------------------------------------
// Dijagnostika 2: jacina krsenja ABC uslova po JEDNOM granicnom trouglu.
// Robin kriterijum (sa e^{+iwt}): n x H + Y E_t = 0,
//     Y = sqrt(epsilon_complex / mu) = k / (omega mu).
// Diskretni FEM Robin je SLAB uslov: pointwise residual ne mora biti nula.
// Pisemo posebne CSV fajlove po rank-u da ne dupliramo elemente.
struct ABCStrongDiagnostics
{
   double residual2[2]={0.,0.};
   double target2[2]={0.,0.};
   double flux_h[2]={0.,0.};
   double flux_e[2]={0.,0.};
   double area[2]={0.,0.};
   long long nfaces[2]={0,0};
};

bool InspectABCFaces(mfem::ParMesh &mesh,
                     const mfem::ParComplexGridFunction &e,
                     double omega,double epsilon0,double mu0,
                     MPI_Comm comm,int rank,bool tem_port,
                     ABCStrongDiagnostics &glob)
{
   const std::string suffix=tem_port?"TEM":"fixed";
   std::ofstream csv("fem_ABC_lokalno_"+suffix+"_rank"+
                     std::to_string(rank)+".csv");
   if(csv)
   {
      csv<<"rank,lokalni_bdr_elem,oznaka,x_mm,y_mm,z_mm,povrsina_mm2,"
            "rms_ABC_rezidual,fluks_H_W,fluks_ABC_E_W,razlika_fluksa_W\n";
      csv<<std::scientific<<std::setprecision(12);
   }
   double local[10]={0.},tot[10]={0.};
   long long local_n[2]={0,0}, global_n[2]={0,0};
   bool ok=true;
   mfem::Vector er(3),ei(3),cr(3),ci(3),hr(3),hi(3),normal(3),xyz(3);
   mfem::Array<int> verts;
   const std::complex<double> eps203(epsilon0*57.55,-1.95/omega);
   const std::complex<double> eps204(epsilon0*2.60,0.);
   const std::complex<double> Y[2]={std::sqrt(eps203/mu0),
                                     std::sqrt(eps204/mu0)};
   for(int b=0;b<mesh.GetNBE();++b)
   {
      const int attr=mesh.GetBdrAttribute(b);
      if(attr!=spolja && attr!=kateter_ulaz)continue;
      const int idx=attr==spolja?0:1;
      mfem::FaceElementTransformations *ft=mesh.GetBdrFaceTransformations(b);
      if(!ft || !ft->Face || !ft->Elem1){ok=false;break;}
      const mfem::IntegrationRule &ir=mfem::IntRules.Get(mesh.GetBdrElementGeometry(b),8);
      double r2=0.,t2=0.,ph=0.,pe=0.,face_area=0.;
      double center[3]={0.,0.,0.};
      int elem=-1,info=-1;
      mesh.GetBdrElementAdjacentElement(b,elem,info);
      if(elem<0 || elem>=mesh.GetNE()){ok=false;break;}
      mesh.GetElementVertices(elem,verts);
      if(verts.Size()!=4){ok=false;break;}
      double ec[3]={0.,0.,0.};
      for(int j=0;j<4;++j)
      {
         const mfem::real_t *x=mesh.GetVertex(verts[j]);
         for(int k=0;k<3;++k)ec[k]+=0.25*x[k];
      }
      for(int j=0;j<ir.GetNPoints();++j)
      {
         const mfem::IntegrationPoint &ip=ir.IntPoint(j);
         ft->Face->SetIntPoint(&ip);
         mfem::IntegrationPoint ipe;
         ft->Loc1.Transform(ip,ipe);
         ft->Elem1->SetIntPoint(&ipe);
         e.real().GetVectorValue(*ft->Elem1,ipe,er);
         e.imag().GetVectorValue(*ft->Elem1,ipe,ei);
         e.real().GetCurl(*ft->Elem1,cr);
         e.imag().GetCurl(*ft->Elem1,ci);
         mfem::CalcOrtho(ft->Face->Jacobian(),normal);
         ft->Face->Transform(ip,xyz);
         double ori=0.;
         for(int k=0;k<3;++k)ori+=normal[k]*(xyz[k]-ec[k]);
         if(!std::isfinite(ori)||std::abs(ori)<1e-24){ok=false;break;}
         if(ori<0.)normal*=-1.;
         const double nn=normal.Norml2(),dA=nn*ip.weight;
         double n[3],etr[3],eti[3],hrr[3],hii[3];
         double en_re=0.,en_im=0.;
         for(int k=0;k<3;++k)
         {
            n[k]=normal[k]/nn;
            en_re+=n[k]*er[k]; en_im+=n[k]*ei[k];
            hrr[k]=-ci[k]/(omega*mu0);
            hii[k]= cr[k]/(omega*mu0);
         }
         for(int k=0;k<3;++k)
         {
            etr[k]=er[k]-en_re*n[k];
            eti[k]=ei[k]-en_im*n[k];
         }
         double residual_point=0.,target_point=0.;
         const double Yr=Y[idx].real(),Yi=Y[idx].imag();
         for(int k=0;k<3;++k)
         {
            const int a=(k+1)%3,c=(k+2)%3;
            const double nxHr=n[a]*hrr[c]-n[c]*hrr[a];
            const double nxHi=n[a]*hii[c]-n[c]*hii[a];
            const double tarRe=Yr*etr[k]-Yi*eti[k];
            const double tarIm=Yr*eti[k]+Yi*etr[k];
            residual_point+=std::pow(nxHr+tarRe,2)+std::pow(nxHi+tarIm,2);
            target_point+=tarRe*tarRe+tarIm*tarIm;
         }
         const double sx=0.5*(er[1]*hrr[2]-er[2]*hrr[1]+ei[1]*hii[2]-ei[2]*hii[1]);
         const double sy=0.5*(er[2]*hrr[0]-er[0]*hrr[2]+ei[2]*hii[0]-ei[0]*hii[2]);
         const double sz=0.5*(er[0]*hrr[1]-er[1]*hrr[0]+ei[0]*hii[1]-ei[1]*hii[0]);
         const double ph_point=(sx*n[0]+sy*n[1]+sz*n[2])*dA;
         double Et2=0.;
         for(int k=0;k<3;++k)Et2+=etr[k]*etr[k]+eti[k]*eti[k];
         const double pe_point=0.5*Yr*Et2*dA;
         r2+=residual_point*dA;
         t2+=target_point*dA;
         ph+=ph_point;pe+=pe_point;
         face_area+=dA;
         for(int k=0;k<3;++k)center[k]+=xyz[k]*dA;
      }
      if(!ok)break;
      if(!(face_area>0.)||!std::isfinite(r2)||!std::isfinite(t2))
      {ok=false;break;}
      local[idx*5]+=r2;
      local[idx*5+1]+=t2;
      local[idx*5+2]+=ph;
      local[idx*5+3]+=pe;
      local[idx*5+4]+=face_area;
      local_n[idx]++;
      if(csv)
      {
         csv<<rank<<','<<b<<','<<attr<<','
            <<1e3*center[0]/face_area<<','
            <<1e3*center[1]/face_area<<','
            <<1e3*center[2]/face_area<<','
            <<1e6*face_area<<','
            <<std::sqrt(r2/std::max(t2,1e-300))<<','
            <<ph<<','<<pe<<','<<ph-pe<<'\n';
      }
   }
   if(!ParallelAll(comm,ok))return false;
   MPI_Allreduce(local,tot,10,MPI_DOUBLE,MPI_SUM,comm);
   MPI_Allreduce(local_n,global_n,2,MPI_LONG_LONG,MPI_SUM,comm);
   for(int j=0;j<2;++j)
   {
      glob.residual2[j]=tot[j*5];glob.target2[j]=tot[j*5+1];
      glob.flux_h[j]=tot[j*5+2];glob.flux_e[j]=tot[j*5+3];
      glob.area[j]=tot[j*5+4];glob.nfaces[j]=global_n[j];
   }
   return true;
}

void WriteABCSummary(const ABCStrongDiagnostics &g,int rank,bool tem_port)
{
   if(rank!=0)return;
   const std::string suf=tem_port?"TEM":"fixed";
   std::ofstream out("fem_ABC_pregled_"+suf+".csv");
   if(out)out<<"oznaka,broj_trouglova,ukupna_povrsina_m2,"
                "L2_relativni_rezidual,P_H_W,P_ABC_E_W,razlika_W\n";
   std::cout<<"\n=== JACINA ODBOJNOG ABC USLOVA PO GRANICI ===\n";
   for(int j=0;j<2;++j)
   {
      const double rel=std::sqrt(g.residual2[j]/std::max(g.target2[j],1e-300));
      std::cout<<"  "<<(j==0?203:204)<<" : RMS |n x H + Y Et| / |Y Et| = "
               <<rel<<", P_H="<<g.flux_h[j]<<" W, P_ABC(E)="
               <<g.flux_e[j]<<" W, ntr="<<g.nfaces[j]<<'\n';
      if(out)out<<(j==0?203:204)<<','<<g.nfaces[j]<<','<<g.area[j]<<','
                 <<rel<<','<<g.flux_h[j]<<','<<g.flux_e[j]<<','
                 <<g.flux_h[j]-g.flux_e[j]<<'\n';
   }
   std::cout<<"  Po-trouglovne vrednosti su u fem_ABC_lokalno_"<<suf
            <<"_rank*.csv (sortirati po rms_ABC_rezidual ili |razlika_fluksa_W|).\n"
            <<"  NAPOMENA: Robin je SLAB uslov, stoga lokalni strong residual NIJE "
               "istovremeno greska solvera.\n"
            <<"================================================\n";
}

// Dijagnostika 3: poredjenje analiticke i diskretne TEM norme, kao i
// procenat poprecnog polja koji ostaje posle izdvajanja TEM moda.
void WritePortModalDetail(const BoundaryPowerFlux &f,double epsilon0,double mu0,
                          int rank,bool tem_port)
{
   if(rank!=0)return;
   const double pi=std::acos(-1.);
   const double n_exact=2.*pi/std::log(coax_r_outer/coax_r_inner);
   const double area_exact=pi*(coax_r_outer*coax_r_outer-coax_r_inner*coax_r_inner);
   const double pinc_nominal=excitation_voltage*excitation_voltage/(2.*f.z0);
   const double E_modal=std::norm(f.ve)*f.mode_norm;
   const double H_modal=std::norm(f.vh)*f.mode_norm;
   const double E_perp=std::max(0.,f.port_et_l2-E_modal);
   const double H_perp=std::max(0.,f.port_eta_ht_l2-H_modal);
   const double e_frac=E_perp/std::max(f.port_et_l2,1e-300);
   const double h_frac=H_perp/std::max(f.port_eta_ht_l2,1e-300);
   std::cout<<"\n=== PRECIZNOST TEM PROJEKCIJE PORTA 201 ===\n"
            <<"  port trouglovi: "<<f.faces[0]<<'\n'
            <<"  povrsina FEM [m^2]: "<<f.area[0]
            <<" ; analiticki anulus [m^2]: "<<area_exact<<'\n'
            <<"  ||TEM(1V)||^2 FEM: "<<f.mode_norm
            <<" ; analiticki: "<<n_exact<<'\n'
            <<"  relat. razlika norme [%]: "<<100.*(f.mode_norm/n_exact-1.)<<'\n'
            <<"  L2 ostatak poprecnog E van osnovnog TEM [%]: "<<100.*e_frac<<'\n'
            <<"  L2 ostatak eta*H van osnovnog TEM [%]: "<<100.*h_frac<<'\n'
            <<"  (L2 procenti su kvadratne norme, nisu direktno procenti snage.)\n"
            <<"  zadati 1 Vpeak => analiticki P_inc [W]: "<<pinc_nominal<<'\n'
            <<"  P_inc iz projektovanog polja [W]: "<<f.modal_inc<<'\n'
            <<"  |V_inc-1 V| [V]: "<<std::abs(f.vinc-std::complex<double>(excitation_voltage,0.))<<'\n'
            <<"  modalni neto [W]: "<<f.modal_net
            <<" ; Poynting neto [W]: "<<-f.out[0]<<'\n'
            <<"  OVO NIJE studija konvergencije PORT mreze: potrebno je ponoviti "
              "racun na guscoj mrezi porta.\n"
            <<"================================================\n";
   std::ofstream out("fem_TEM_projekcija_provera_"+
                     std::string(tem_port?"TEM":"fixed")+".csv");
   if(out)
   {
      out<<std::scientific<<std::setprecision(13)
         <<"rezim,N_fem,N_analiticki,area_fem_m2,area_analiticka_m2,"
           "E_L2_nemodalno_frakcija,H_L2_nemodalno_frakcija,"
           "Vinc_re,Vinc_im,Pinc_zadata_1V_W,Pinc_projekcija_W,"
           "Pnet_modal_W,Pnet_poynting_W,S11_abs\n"
         <<(tem_port?"TEM":"fixed")<<','<<f.mode_norm<<','<<n_exact<<','
         <<f.area[0]<<','<<area_exact<<','<<e_frac<<','<<h_frac<<','
         <<f.vinc.real()<<','<<f.vinc.imag()<<','<<pinc_nominal<<','
         <<f.modal_inc<<','<<f.modal_net<<','<<-f.out[0]<<','
         <<std::abs(f.s11)<<'\n';
   }
}


} // namespace

int main(int argc, char *argv[])
{
   mfem::Mpi::Init(argc, argv);
   const MPI_Comm comm = MPI_COMM_WORLD;
   const int rank = mfem::Mpi::WorldRank();
   const int np = mfem::Mpi::WorldSize();
   mfem::Hypre::Init();
   if (rank == 0)
   {
      std::cerr << "*** MPI MUMPS ANTENA + TUMOR 107 (DIREKTNI RESAVAC) ***" << std::endl;
   }

   // --port fixed: istorijski 1 V Dirichlet (referentni rezultat)
   // --port tem: incidentni 1 Vpeak TEM + izlazni Robin uslov.
   if (argc < 2)
   {
      if (rank == 0)
      {
         std::cerr << "Upotreba: mpirun -np 2 ./antena_em_p putanja.msh "
                   << "[--rtol TOL] [--port fixed|tem]\n";
      }
      return 1;
   }
   double true_tolerance=1e-8;
   bool tem_port=false;
   if((argc-2)%2!=0)
   {
      if(rank==0)std::cerr<<"Opcije: [--rtol TOL] [--port fixed|tem]\n";
      return 1;
   }
   for(int k=2;k<argc;k+=2)
   {
      const std::string key=argv[k],value=argv[k+1];
      if(key=="--port")
      {
         if(value=="tem")tem_port=true;
         else if(value=="fixed")tem_port=false;
         else {if(rank==0)std::cerr<<"Port mora biti fixed ili tem.\n";return 1;}
      }
      else if(key=="--rtol")
      {
         try {true_tolerance=std::stod(value);}
         catch(const std::exception &){if(rank==0)std::cerr<<"Neispravan rtol.\n";return 1;}
      }
      else {if(rank==0)std::cerr<<"Nepoznata opcija: "<<key<<'\n';return 1;}
   }
   if (!std::isfinite(true_tolerance) || true_tolerance <= 0.0 ||
       true_tolerance >= 1.0)
   {
      if(rank==0)std::cerr<<"Neispravna tolerancija.\n";
      return 1;
   }

   std::unique_ptr<mfem::ParMesh> pmesh;
   {
      std::ifstream input(argv[1], std::ios::binary);
      if (!ParallelAll(comm, static_cast<bool>(input)))
      {
         if (rank == 0) { std::cerr << "Ne mogu da otvorim mrezu: " << argv[1] << '\n'; }
         return 1;
      }
      mfem::Mesh cela_mreza(input, 1, 0);
      if (cela_mreza.Dimension() != 3 || cela_mreza.SpaceDimension() != 3)
      {
         if (rank == 0) { std::cerr << "Ocekujem 3D volumensku mrezu.\n"; }
         return 2;
      }
      for (int i = 0; i < cela_mreza.attributes.Size(); ++i)
      {
         const int attr = cela_mreza.attributes[i];
         if (attr < metal || attr > tumor)
         {
            if (rank == 0)
            {
               std::cerr << "Nepoznata zapremina " << attr
                         << " (dozvoljene oznake su 101 do 107).\n";
            }
            return 2;
         }
      }
      mfem::Array<int> em_atributi;
      for (int attr : {izolator, slot1, slot2, kateter, tkivo, tumor})
      {
         if (!HasAttribute(cela_mreza.attributes, attr))
         {
            if (rank == 0) { std::cerr << "Nedostaje zapremina " << attr << '\n'; }
            return 2;
         }
         em_atributi.Append(attr);
      }

      mfem::SubMesh mesh = mfem::SubMesh::CreateFromDomain(cela_mreza, em_atributi);
      const bool valid = ValidateSerialMesh(cela_mreza, mesh, rank);
      if (!ParallelAll(comm, valid)) { return 2; }
      if (rank == 0)
      {
         std::cout << "Pravim paralelnu mrezu za " << np << " MPI procesa...\n";
      }
      pmesh.reset(new mfem::ParMesh(comm, mesh));
   } // serijska mreza vise nije potrebna

   mfem::ParMesh &mesh = *pmesh;
   const int local_elements = mesh.GetNE();
   int global_elements = 0;
   MPI_Allreduce(&local_elements, &global_elements, 1, MPI_INT, MPI_SUM, comm);
   if (rank == 0)
   {
      std::cout << "Globalni broj EM elemenata: " << global_elements << '\n';
   }

   // MPI dijagnostika: konstrukcija ND prostora i globalnih DOF.
   std::cerr << "[rank " << rank << "] PRE ND_FECollection" << std::endl;
   mfem::ND_FECollection fec(1, mesh.Dimension());
   std::cerr << "[rank " << rank << "] PRE ParFiniteElementSpace" << std::endl;
   mfem::ParFiniteElementSpace fes(&mesh, &fec);
   std::cerr << "[rank " << rank << "] POSLE ParFiniteElementSpace / PRE GlobalTrueVSize" << std::endl;

   // GlobalTrueVSize moze zahtevati komunikaciju svih MPI procesa.
   const HYPRE_BigInt global_nd_dofs = fes.GlobalTrueVSize();
   std::cerr << "[rank " << rank << "] POSLE GlobalTrueVSize = "
             << global_nd_dofs << std::endl;
   if (rank == 0)
   {
      std::cout << "Globalni broj kompleksnih ND DOF: "
                << global_nd_dofs << '\n'
                << "Globalna dimenzija realnog blok-sistema: "
                << 2 * global_nd_dofs << std::endl;
   }

   // Neki MPI proces moze nemati nijedan element odredjene oznake.
   // Markeri ipak moraju biti dovoljni za SVE globalne ID-jeve.
   const int max_bdr_attr = std::max(kateter_ulaz,
       mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0);
   mfem::Array<int> port_marker(max_bdr_attr), pec_marker(max_bdr_attr);
   mfem::Array<int> ess_marker(max_bdr_attr);
   port_marker = 0;
   pec_marker = 0;
   ess_marker = 0;
   port_marker[port - 1] = 1;
   pec_marker[pec - 1] = 1;
   ess_marker[port - 1] = tem_port ? 0 : 1;
   ess_marker[pec - 1] = 1;
   mfem::Array<int> ess_tdof;
   fes.GetEssentialTrueDofs(ess_marker, ess_tdof);

   mfem::ParComplexGridFunction e(&fes);
   e = std::complex<mfem::real_t>(0.0, 0.0);
   mfem::VectorFunctionCoefficient port_e(3, PortE);
   mfem::Vector zero_vec(3);
   zero_vec = 0.0;
   mfem::VectorConstantCoefficient zero(zero_vec);
   if(!tem_port)e.ProjectBdrCoefficientTangent(port_e, zero, port_marker);
   e.ProjectBdrCoefficientTangent(zero, zero, pec_marker);

   const double pi = std::acos(-1.0);
   const double frequency = 2.45e9;
   const double omega = 2.0 * pi * frequency;
   const double epsilon0 = 8.8541878128e-12;
   const double mu0 = 4.0 * pi * 1e-7;

   struct Material
   {
      int attr;
      const char *name;
      double eps_r;
      double sigma;
      double mu_r;
   };
   const Material materials[] = {
      {izolator, "Izolator", 2.03, 0.0, 1.0},
      {slot1, "Slot1", 1.00, 0.0, 1.0},
      {slot2, "Slot2", 1.00, 0.0, 1.0},
      {kateter, "Kateter", 2.60, 0.0, 1.0},
      {tkivo, "Tkivo", 57.55, 1.95, 1.0},
      // O'Rourke et al., Phys Med Biol 52 (2007), tabela 3,
      // in vivo maligno tkivo jetre na 2.45 GHz, srednje vrednosti.
      {tumor, "Tumor", 62.44, 2.18, 1.0}
   };

   const int nattr = std::max(tumor,
       mesh.attributes.Size() ? mesh.attributes.Max() : 0);
   mfem::Vector curl_values(nattr), mass_values(nattr), loss_values(nattr);
   curl_values = 0.0;
   mass_values = 0.0;
   loss_values = 0.0;
   if (rank == 0) { std::cout << "Frekvencija: " << frequency/1e9 << " GHz\n"; }
   for (const Material &m : materials)
   {
      const int idx = m.attr - 1;
      curl_values(idx) = 1.0 / (mu0 * m.mu_r);
      mass_values(idx) = -omega * omega * epsilon0 * m.eps_r;
      loss_values(idx) = omega * m.sigma;
      if (rank == 0)
      {
         std::cout << "  " << m.attr << " " << m.name << ": eps_r="
                   << m.eps_r << ", sigma=" << m.sigma << " S/m\n";
      }
   }
   mfem::PWConstCoefficient curl_coef(curl_values);
   mfem::PWConstCoefficient mass_coef(mass_values);
   mfem::PWConstCoefficient loss_coef(loss_values);

   // Isti ABC prvog reda kao u serijskoj verziji (nije PML).
   mfem::Array<int> abc_marker(max_bdr_attr);
   abc_marker = 0;
   abc_marker[spolja - 1] = 1;
   abc_marker[kateter_ulaz - 1] = 1;
   mfem::Vector abc_real_values(max_bdr_attr), abc_imag_values(max_bdr_attr);
   abc_real_values = 0.0;
   abc_imag_values = 0.0;
   const std::complex<double> I(0.0, 1.0);
   for (const Material &m : materials)
   {
      if (m.attr != tkivo && m.attr != kateter) { continue; }
      const int boundary_attr = (m.attr == tkivo) ? spolja : kateter_ulaz;
      const double mu = mu0 * m.mu_r;
      const std::complex<double> epsilon_complex(epsilon0*m.eps_r,
                                                 -m.sigma/omega);
      const std::complex<double> k = omega*std::sqrt(mu*epsilon_complex);
      const std::complex<double> q = I*k/mu;
      abc_real_values(boundary_attr-1) = q.real();
      abc_imag_values(boundary_attr-1) = q.imag();
   }
   const double k_coax=omega*std::sqrt(mu0*epsilon0*coax_epsr);
   if(tem_port)
   {
      // Robin: n x mu^-1 curl E = i k/mu (E_t - 2E_inc,t).
      // Proveriti na izlazu pomocu modalne dekompozicije.
      abc_marker[port-1]=1;
      abc_imag_values(port-1)=k_coax/mu0;
      tem_forcing_k_over_mu=k_coax/mu0;
   }
   if(rank==0)
   {
      std::cout<<"Port 201: "<<(tem_port?"incidentni TEM 1 Vpeak + Robin":"Dirichlet ukupni E 1 Vpeak")<<'\n';
      std::cout<<"Koaksijalni Z0 [ohm]: "
               <<std::sqrt(mu0/(epsilon0*coax_epsr))*std::log(coax_r_outer/coax_r_inner)/(2.*pi)<<'\n';
   }
   mfem::PWConstCoefficient abc_real_coef(abc_real_values);
   mfem::PWConstCoefficient abc_imag_coef(abc_imag_values);

   // Fizicki operator je NEPROMENJEN u odnosu na serijski model.
   mfem::ParSesquilinearForm a(&fes, conv);
   a.AddDomainIntegrator(new mfem::CurlCurlIntegrator(curl_coef), nullptr);
   a.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(mass_coef),
                         new mfem::VectorFEMassIntegrator(loss_coef));
   a.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(abc_real_coef),
                           new mfem::VectorFEMassIntegrator(abc_imag_coef),
                           abc_marker);
   if (rank == 0) { std::cout << "Sastavljam paralelni EM operator...\n"; }
   a.Assemble();
   a.Finalize();
   // Kod ParSesquilinearForm postavlja se politika na sastavnim formama.
   a.real().SetDiagonalPolicy(mfem::Operator::DIAG_ONE);
   a.imag().SetDiagonalPolicy(mfem::Operator::DIAG_ZERO);

   mfem::ParComplexLinearForm b(&fes, conv);
   mfem::VectorFunctionCoefficient forcing_im(3,PortTEMImagForcing);
   if(tem_port)
   {
      // VectorFEBoundaryTangentLFIntegrator racuna (n x f) . v.
      // Imaginarni koeficijent f je iz PortTEMImagForcing.
      // Koeficijent mora ostati ziv DO kraja Assemble().
      b.AddBoundaryIntegrator(nullptr,
          new mfem::VectorFEBoundaryTangentLFIntegrator(forcing_im),port_marker);
      b.Assemble();
   }
   else {b=0.0;}
   mfem::OperatorHandle A;
   mfem::Vector X, B;
   a.FormLinearSystem(ess_tdof, e, b, A, X, B);
   const int n = fes.GetTrueVSize(); // lokalni kompleksni ND DOF

   const double rhs_norm = ParallelNorm(comm, B);
   if (rank == 0)
   {
      std::cout << std::scientific << std::setprecision(8)
                << "Globalna norma ||B||_2 = " << rhs_norm << '\n';
   }
   if (!std::isfinite(rhs_norm) || rhs_norm == 0.0)
   {
      if (rank == 0) { std::cerr << "Neispravna desna strana.\n"; }
      return 4;
   }

   // MUMPS prihvata HypreParMatrix, a nas A je kompleksan operator.
   // MFEM formira realni 2x2 blok-sistem sa istom lokalnom raspodelom
   // [Re(E), Im(E)] i konvencijom HERMITIAN kao originalni operator.
   auto *complex_A = dynamic_cast<mfem::ComplexOperator*>(A.Ptr());
   if (!ParallelAll(comm, complex_A != nullptr))
   {
      if (rank == 0) { std::cerr << "Operator A nije ComplexOperator.\n"; }
      return 5;
   }
   std::unique_ptr<mfem::ComplexHypreParMatrix> complex_hypre(
       complex_A->AsComplexHypreParMatrix());
   std::unique_ptr<mfem::HypreParMatrix> real_system(
       complex_hypre->GetSystemMatrix());
   if (!ParallelAll(comm, real_system != nullptr))
   {
      if (rank == 0) { std::cerr << "Ne mogu da formiram MUMPS matricu.\n"; }
      return 5;
   }

   // Nezavisna provera da je MUMPS matrica IDENTICNA originalnom
   // operatoru i da je raspored realnog i imaginarnog dela kompatibilan.
   mfem::Vector probe(B.Size()), y_original(B.Size()), y_matrix(B.Size());
   for (int i = 0; i < probe.Size(); ++i)
   {
      probe[i] = std::sin(0.011 * (i + 1) * (rank + 1));
   }
   A->Mult(probe, y_original);
   real_system->Mult(probe, y_matrix);
   y_matrix -= y_original;
   const double operator_diff = ParallelNorm(comm, y_matrix) /
       std::max(1.0, ParallelNorm(comm, y_original));
   if (rank == 0)
   {
      std::cout << "Provera konverzije kompleksne matrice: "
                << operator_diff << '\n';
   }
   if (!std::isfinite(operator_diff) || operator_diff > 1e-10)
   {
      if (rank == 0)
      {
         std::cerr << "MUMPS matrica se razlikuje od originalnog A.\n";
      }
      return 5;
   }

   mfem::Vector residual(B.Size());
   A->Mult(X, residual);
   residual -= B;
   const double initial_norm = ParallelNorm(comm, residual) / rhs_norm;
   if (rank == 0)
   {
      std::cout << "Pocetni stvarni relativni ostatak: " << initial_norm << '\n'
                << "Pocinjem MPI MUMPS: nproc=" << np << '\n';
   }

   // Nema HypreAMS preduslovljavaca ni GMRES iteracija: MUMPS pravi
   // simbolicku i numericku LU faktorizaciju distribuirane matrice.
   // Matriza nije pozitivno definitna: siguran izbor je UNSYMMETRIC.
   mfem::MUMPSSolver mumps(comm);
   mumps.SetPrintLevel(2);
   mumps.SetMatrixSymType(mfem::MUMPSSolver::UNSYMMETRIC);
   mumps.SetReorderingStrategy(mfem::MUMPSSolver::AUTOMATIC);
   mumps.SetOperator(*real_system); // analiza + faktorizacija, kolektivno
   mumps.Mult(B, X);                // resavanje, kolektivno

   // Nezavisna MPI provera ORIGINALNOG sistema (ne preduslovljenog).
   A->Mult(X, residual);
   residual -= B;
   const double relative_residual = ParallelNorm(comm, residual) / rhs_norm;
   bool local_finite = true;
   for (int i = 0; i < X.Size(); ++i)
   {
      if (!std::isfinite(X[i])) { local_finite = false; break; }
   }
   const bool finite_solution = ParallelAll(comm, local_finite);

   double local_boundary_error = 0.0;
   double local_boundary_scale = 1.0;
   for (int i = 0; i < ess_tdof.Size(); ++i)
   {
      const int dof = ess_tdof[i];
      for (int block = 0; block < 2; ++block)
      {
         const int k = dof + block*n;
         local_boundary_error = std::max(local_boundary_error,
                                         static_cast<double>(std::abs(X[k]-B[k])));
         local_boundary_scale = std::max(local_boundary_scale,
                                         static_cast<double>(std::abs(B[k])));
      }
   }
   const double boundary_error = ParallelMax(comm, local_boundary_error);
   const double boundary_scale = ParallelMax(comm, local_boundary_scale);
   const bool accepted = finite_solution &&
                         std::isfinite(relative_residual) &&
                         relative_residual <= true_tolerance &&
                         boundary_error <= 1e-10*boundary_scale;
   if (rank == 0)
   {
      std::cout << std::scientific << std::setprecision(10)
                << "\nMPI MUMPS rezultat\n"
                << "  Stvarni ||AX-B||/||B||: " << relative_residual << '\n'
                << "  Zahtevano: " << true_tolerance << '\n'
                << "  Max odstupanje zadatih DOF: " << boundary_error << '\n'
                << "  Polje prihvaceno: " << (accepted?"DA":"NE") << '\n';
   }
   if (!accepted)
   {
      if (rank == 0)
      {
         std::cerr << "Resavac nije prosao proveru; novo polje nije sacuvano.\n";
      }
      return 5;
   }

   a.RecoverFEMSolution(X, b, e);

   // sigma se cita iz ISTE tablice koja se koristi u Maxwellovom operatoru,
   // kako se numericka integracija ne bi razlikovala po materijalima.
   double sigma106 = -1.0, sigma107 = -1.0;
   for (const Material &m : materials)
   {
      if (m.attr == tkivo) { sigma106 = m.sigma; }
      if (m.attr == tumor) { sigma107 = m.sigma; }
   }
   if (sigma106 < 0.0 || sigma107 < 0.0)
   {
      if (rank == 0) { std::cerr << "Nedostaje sigma za 106/107.\n"; }
      return 6;
   }

   PowerIntegrals integral;
   if (!IntegrateTissuePower(mesh, e, sigma106, sigma107, comm, integral))
   {
      if (rank == 0)
      {
         std::cerr << "Greska u MFEM zapreminskoj integraciji Q.\n";
      }
      return 6;
   }
   PrintPower(integral, rank, tem_port);
   if (rank == 0)
   {
      if (SavePowerCSV(integral, sigma106, sigma107, tem_port))
      {
         std::cout << "Sacuvan CSV sa FEM integralima, rezim " << (tem_port?"TEM":"fixed") << "\n";
      }
      else
      {
         std::cerr << "UPOZORENJE: ne mogu da sacuvam CSV integrala.\n";
      }
   }
   // Za isti reseni FEM sistem proveravamo razlicite kvadrature i
   // obradu boundary vs volume H(curl) evaluacije.
   const std::array<int,4> orders={{2,4,6,8}};
   BoundaryPowerFlux reference;
   std::ofstream diag;
   if(rank==0)
   {
      diag.open(tem_port?"fem_kvadratura_H_TEM.csv":"fem_kvadratura_H_fixed.csv");
      if(diag)diag<<"metod,kvadratura,port201_W,pec202_W,abc203_W,abc204_W,"
                     "rezidual_W,rezidual_proc,abc203_iz_E_W,abc204_iz_E_W,"
                     "P_TEM_inc_W,P_TEM_ref_W,P_TEM_net_W,S11_abs\n";
      std::cout<<"\n=== KVADRATURA I EVALUACIJA H = i curl(E)/(omega mu) ===\n";
   }
   for(bool volume_side : {false,true})
   {
      for(int qorder:orders)
      {
         BoundaryPowerFlux flux;
         if(!IntegratePoyntingFlux(mesh,e,omega,mu0,epsilon0,comm,
                                   qorder,volume_side,flux))
         {
            if(rank==0)std::cerr<<"Greska granicne evaluacije za kvadraturu "
                                  <<qorder<<".\n";
            return 7;
         }
         const PowerBalance pb=CalculatePowerBalance(integral,flux);
         if(rank==0)
         {
            std::cout<<"  "<<(volume_side?"volume_side":"boundary")
                     <<"  Q="<<qorder
                     <<" P201net [W]="<<pb.pin_net
                     <<"  rezidual [W]="<<pb.residual
                     <<"  odstupanje [%]="<<100.*pb.relative<<'\n';
            if(diag)diag<<(volume_side?"volume_side":"boundary")<<','
               <<qorder<<','<<flux.out[0]<<','<<flux.out[1]<<','
               <<flux.out[2]<<','<<flux.out[3]<<','<<pb.residual<<','
               <<100.*pb.relative<<','<<flux.abc_e203<<','<<flux.abc_e204
               <<','<<flux.modal_inc<<','<<flux.modal_ref<<','<<flux.modal_net
               <<','<<std::abs(flux.s11)<<'\n';
         }
         if(volume_side && qorder==8)reference=flux;
      }
   }
   const PowerBalance balance=CalculatePowerBalance(integral,reference);
   PrintPowerBalance(integral,reference,balance,rank,tem_port);

   WeakEnergyDiagnostics weak;
   if(!CalculateWeakEnergy(*A.Ptr(),X,B,omega,epsilon0,mu0,
                           reference,integral,comm,tem_port,weak))
   {
      if(rank==0)std::cerr<<"Greska slabe energetske provere.\n";
      return 8;
   }
   WriteWeakEnergy(weak,reference,rank,tem_port);

   ABCStrongDiagnostics strong_abc;
   if(!InspectABCFaces(mesh,e,omega,epsilon0,mu0,comm,rank,tem_port,strong_abc))
   {
      if(rank==0)std::cerr<<"Greska lokalne ABC dijagnostike.\n";
      return 8;
   }
   WriteABCSummary(strong_abc,rank,tem_port);
   WritePortModalDetail(reference,epsilon0,mu0,rank,tem_port);
   if(rank==0)
   {
      const double s11abs=std::abs(reference.s11);
      const double s11db=s11abs>0?20.*std::log10(s11abs): -std::numeric_limits<double>::infinity();
      std::cout<<"\n=== KOAKSIJALNI TEM MODALNI PORT 201 ===\n"
               <<"Rezim: "<<(tem_port?"TEM Robin sa 1 Vpeak incidentnog talasa"
                                  :"Dirichlet ukupno E: modovi su dijagnostika")<<'\n'
               <<"Z0 [ohm]: "<<reference.z0<<'\n'
               <<"V_E (kompleksno) [V]: "<<reference.ve<<'\n'
               <<"V_H (kompleksno) [V]: "<<reference.vh<<'\n'
               <<"V_inc [V]: "<<reference.vinc<<'\n'
               <<"V_ref [V]: "<<reference.vref<<'\n'
               <<"P_inc TEM [W]: "<<reference.modal_inc<<'\n'
               <<"P_ref TEM [W]: "<<reference.modal_ref<<'\n'
               <<"P_inc - P_ref [W]: "<<reference.modal_net<<'\n'
               <<"P_neto flux 201 [W]: "<<balance.pin_net<<'\n'
               <<"Razlika TEM-neto - flux-port [W]: "
               <<reference.modal_net-balance.pin_net<<'\n'
               <<"S11 = V_ref/V_inc: "<<reference.s11<<'\n'
               <<"|S11|: "<<s11abs<<"; S11 [dB]: "<<s11db<<'\n'
               <<"Reflektovana frakcija TEM: "<<s11abs*s11abs<<'\n';
      const double modal_mismatch=std::abs(reference.modal_net-balance.pin_net)/
         std::max(std::abs(balance.pin_net),1e-30);
      if(modal_mismatch>0.05)
         std::cout<<"UPOZORENJE: TEM neto i ukupni Poynting port se razlikuju >5%. "
                   "S11 interpretirati OPREZNO (visoki/evanescentni modovi, "
                   "diskretizacija, polozaj porta).\n";
      // Skaliranje na 10 W nominalne incidentne pobude se radi
      // kasnije iz zadatog TEM napona i Z0 (NE iz modalne projekcije).
      if(tem_port)
         std::cout<<"NOMINALNA 10W pobuda: faktor ce biti izracunat iz "
                     "V_inc_peak^2/(2 Z0), ne iz P_inc projekcije.\n";
      if(tem_port)
      {
         std::cout<<"Port je Robin pobudjen zadatim incidentnim profilom 1 Vpeak; "
                    "izracunati V_inc je projekcija numerickog resenja.\n";
      }
      std::cout<<"Ovo je TEM projekcija, vazeca ako su visi modovi "
                "i evanescentni doprinosi na portu zanemarljivi.\n"
                <<"====================================================\n";
      if(SavePowerBalanceCSV(integral,reference,balance,tem_port))
         std::cout<<"Sacuvan CSV energetskog bilansa: " << (tem_port?"TEM":"fixed") << "\n";
      std::ofstream modal(tem_port?"fem_port_TEM_S11.csv":"fem_port_fixed_S11_dijagnostika.csv");
      if(modal)
      {
         modal<<std::scientific<<std::setprecision(12)
              <<"rezim,Z0_ohm,P_inc_W,P_ref_W,P_net_TEM_W,P_net_flux201_W,"
                "V_inc_re,V_inc_im,V_ref_re,V_ref_im,S11_re,S11_im,S11_abs,S11_dB\n"
              <<(tem_port?"tem":"fixed")<<','<<reference.z0<<','
              <<reference.modal_inc<<','<<reference.modal_ref<<','
              <<reference.modal_net<<','<<balance.pin_net<<','
              <<reference.vinc.real()<<','<<reference.vinc.imag()<<','
              <<reference.vref.real()<<','<<reference.vref.imag()<<','
              <<reference.s11.real()<<','<<reference.s11.imag()<<','
              <<s11abs<<','<<s11db<<'\n';
      }
   }
   // NOMINALNA SNAGA. Ovaj faktor ne zavisi od modalne PROJEKCIJE
   // rekonstruisanog polja, nego od 1 Vpeak i analitickog Z0 na portu 201.
   const double z0_nominal=std::sqrt(mu0/(epsilon0*coax_epsr)) *
      std::log(coax_r_outer/coax_r_inner)/(2.0*pi);
   const double p_inc_nominal=excitation_voltage*excitation_voltage/
      (2.0*z0_nominal);
   const double e_scale=tem_port ?
      std::sqrt(target_incident_power_W/p_inc_nominal) : 1.0;
   const double q_scale=e_scale*e_scale;
   if (!(p_inc_nominal>0.0) || !std::isfinite(e_scale) ||
       !std::isfinite(q_scale))
   {
      if(rank==0)std::cerr<<"Neispravna nominalna TEM normalizacija.\n";
      return 9;
   }

   // L2(0) = po jedna nezavisna scalar DOF po tetraedru (cell-average Q).
   mfem::L2_FECollection q_collection(0,mesh.Dimension());
   mfem::ParFiniteElementSpace q_fes(&mesh,&q_collection);
   mfem::ParGridFunction q_ref(&q_fes);
   std::unique_ptr<mfem::ParGridFunction> q_10W;
   std::unique_ptr<mfem::ParGridFunction> er_10W, ei_10W;
   if(tem_port)
   {
      q_10W.reset(new mfem::ParGridFunction(&q_fes));
      er_10W.reset(new mfem::ParGridFunction(&fes));
      ei_10W.reset(new mfem::ParGridFunction(&fes));
      for (int i=0;i<fes.GetVSize();++i)
      {
         (*er_10W)(i)=e_scale*e.real()(i);
         (*ei_10W)(i)=e_scale*e.imag()(i);
      }
   }
   ThermalSourceVerification srccheck;
   if (!MakeThermalSourceDG0(mesh,e,q_fes,q_ref,q_10W.get(),
                             sigma106,sigma107,q_scale,comm,srccheck))
   {
      if(rank==0)std::cerr<<"Neuspeo DG0 izvoz izvora toplote.\n";
      return 9;
   }
   const double maxdiff=std::max(std::abs(srccheck.p106_ref-integral.p106),
                                 std::abs(srccheck.p107_ref-integral.p107));
   if(maxdiff>1e-8*std::max(1.0,integral.p106+integral.p107))
   {
      if(rank==0)std::cerr<<"Greska: DG0 integral Q nije konzervativan!\n";
      return 9;
   }
   if(rank==0)
   {
      if(tem_port)
      {
         std::cout<<"\n=== IZVOR TOPLOTE: 10 W NOMINALNE INCIDENTNE TEM SNAGE ===\n"
                  <<"  V_inc_peak referentno [V]: "<<excitation_voltage<<'\n'
                  <<"  Z0 [ohm]: "<<z0_nominal<<'\n'
                  <<"  P_inc nominalno referentno [W]: "<<p_inc_nominal<<'\n'
                  <<"  Ciljna incidentna snaga [W]: "<<target_incident_power_W<<'\n'
                  <<"  Skala E [1]: "<<e_scale<<'\n'
                  <<"  Skala Q [1]: "<<q_scale<<'\n'
                  <<"  P_abs 106 pri P_inc=10W [W]: "<<srccheck.p106_10W<<'\n'
                  <<"  P_abs 107 pri P_inc=10W [W]: "<<srccheck.p107_10W<<'\n'
                  <<"  Ukupni apsorbovani W: "<<(srccheck.p106_10W+srccheck.p107_10W)<<'\n'
                  <<"  Q je definisan SAMO u 106/107, ostali EM materijali Q=0.\n"
                  <<"  Polje Q_Pinc_10W_W_m3 je zapreminski DG0 prosek,\n"
                  <<"  ne originalna tackasta raspodela na kvadraturnim tackama.\n"
                  <<"  Nominalno 10W nije isto sto i 10W neto na portu.\n"
                  <<"========================================================\n";
         if(!WriteNominalPower10WCSV(integral,srccheck,z0_nominal,
                                     p_inc_nominal,e_scale,q_scale,tem_port))
            std::cerr<<"UPOZORENJE: CSV sa izvorom toplote nije sacuvan.\n";
      }
      else
      {
         std::cout<<"Q za 1 V ukupni Dirichlet je izvezen, ali bez \n"
                    "10 W incidentnog skaliranja (samo TEM rezim).\n";
      }
   }
   return SaveField(mesh,e,q_ref,q_10W.get(),er_10W.get(),ei_10W.get(),
                    comm,rank,tem_port);
}
