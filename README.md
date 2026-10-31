# Microwave Ablation Research Internship

**Three-Dimensional Finite-Element Modeling of a Patient-Derived Liver Tumor Using a Dual-Slot Coaxial Antenna**

Research internship at the **Institute of Physics Belgrade**, October 2026.  
**Author:** Viktoria Đorđević · **Research mentor:** Dr Nikola Bošković

[Read the research manuscript (PDF)](report/Microwave%20Ablation%20%281%29.pdf)

## Overview

This project investigates electromagnetic energy deposition, transient tissue heating, and model-based thermal injury around a dual-slot coaxial microwave applicator. The workflow combines a three-dimensional Gmsh/OpenCASCADE geometry, an MFEM electromagnetic solver at 2.45 GHz, and a parallel Pennes bioheat solver with Arrhenius injury modeling.

The manuscript also examines kinetic sensitivity and fixed-placement power–duration protocols. The tumor geometry derives from case **3Dircadb1.2** of the **3D-IRCADb-01** dataset, as documented and referenced in the manuscript. The enclosing healthy-tissue cylinder is a computational surrogate, not a segmented whole liver.

## Repository contents

| Path | Contents |
| --- | --- |
| `report/Microwave Ablation (1).pdf` | 31-page internship manuscript, figures, results, and references |
| `src/antena_em_p_TEM_S11_Q_Pinc10W.cpp` | MPI/MUMPS electromagnetic solver, TEM feed, S11 diagnostics, and nominal 10 W heat-source normalization |
| `src/bioheat_R35_unified.cpp` | Unified transient Pennes/Arrhenius solver with power, duration, perfusion, and kinetic parameters |
| `geometry/antena_V3_2_tumor107_R35_port_h010.geo` | Antenna/tumor geometry and local mesh sizing, with 35 mm tissue-domain radius |
| `geometry/tumor_h1_5_107_izvan_katetera_m.brep` | Tumor CAD solid in metres, with catheter envelope removed |
| `geometry/antena_tumor_R35_port_h010_netgen_ascii.msh` | Supplied Gmsh 4.1 ASCII mesh |
| `SHA256SUMS` | SHA-256 checksums for the six unchanged research files |

Keep the `.geo` and `.brep` together: the geometry imports the BREP by its original filename. The BREP already uses metres; the script sets `Geometry.OCCScaling = 1`.

## Computational workflow

1. Construct conforming antenna, catheter, tumor, and surrounding-tissue regions in Gmsh.
2. Solve the complex Maxwell problem with first-order Nédélec elements and MPI/MUMPS.
3. Normalize the TEM solution to 10 W nominal incident power and export a conservative DG0 heat-source field.
4. Transfer cellwise source values to tissue regions 106 and 107 using centroid and volume checks.
5. Solve transient Pennes heat transfer and integrate Arrhenius injury; inspect energy balances, temperature, and model-injury coverage.

Physical volume tags are 101 (metal), 102 (coaxial dielectric), 103/104 (slots), 105 (catheter), 106 (healthy tissue), and 107 (tumor). Boundary tags are 201 (port), 202 (PEC), and 203 (outer absorbing boundary); the EM code identifies boundary 204 on the catheter inlet after submesh extraction.

## Selected manuscript results

These are results reported in the supplied manuscript, **not independently rerun for this repository**.

| Quantity | Reported result |
| --- | --- |
| Frequency | 2.45 GHz |
| Full mesh | 294,380 tetrahedra |
| Electromagnetic submesh | 226,153 tetrahedra |
| Projected port S11 | approximately −11.90 dB |
| Absorbed tissue power at 10 W incident | approximately 9.00497 W |
| Maximum tissue temperature at 10 W, 300 s | 85.728 °C |
| Tumor volume meeting the reference injury criterion at 10 W, 300 s | 30.651% |

The injury criterion is a model-derived fraction of at least 0.98. Its interpretation depends strongly on kinetic parameters. Spatial/time convergence and experimental validation remain outstanding; these numerical outputs are not clinically validated treatment predictions.

## Requirements and execution

The manuscript records **MFEM 4.10** and **MUMPS 5.8.2**, with a C++17 compiler, MPI, and HYPRE. Use compatible MPI installations across dependencies. Gmsh is needed to inspect or regenerate the geometry; ParaView is used to visualize parallel VTK collections. MUMPS factorization can require several GB of memory.

This upload preserves the supplied source files. It does **not** include the original CMake configuration or a tested build environment. Compile and link the two sources against an MFEM installation with MPI/HYPRE/MUMPS enabled, naming the resulting executables `antena_em_p` and `bioheat_R35_unified` for the examples below.

From the repository root, after placing compiled executables in `build/`:

```bash
PROJECT_ROOT="$PWD"
mkdir -p runs/em-10w
cd runs/em-10w
mpirun -np 2 "$PROJECT_ROOT/build/antena_em_p" \
  "$PROJECT_ROOT/geometry/antena_tumor_R35_port_h010_netgen_ascii.msh" \
  --port tem --rtol 1e-8
```

Explicitly pass `--port tem`: the source defaults to the historical fixed-field port. The 10 W incident normalization is defined in the EM source, which has no `--power` option. Run each experiment in a fresh directory because output names are reused.

### Heat-source transfer and thermal execution

The EM solver exports the field `Q_Pinc_10W_W_m3` in a ParaView collection. The diagnostic `fem_izvor_toplote_Pinc10W.csv` is a key–value power summary and **cannot be used as the thermal cell-source input**.

The thermal solver requires a separate `Q_Pinc10W_cells.csv` with columns:

```text
cx_m,cy_m,cz_m,volume_m3,q_W_m3
```

The `export_q_from_paraview.py` script referenced by the manuscript and the precomputed cell-source CSV were **not among the six supplied files**. Obtain that exporter or a validated equivalent before attempting the thermal step. Preserve the correspondence with the exact mesh and check the original integrated source power (default reference: 9.0049672684 W).

Once a validated cell-source CSV is available, an example thermal invocation is:

```bash
mkdir -p "$PROJECT_ROOT/runs/thermal-10w-300s"
cd "$PROJECT_ROOT/runs/thermal-10w-300s"
mpirun -np 2 "$PROJECT_ROOT/build/bioheat_R35_unified" \
  "$PROJECT_ROOT/geometry/antena_tumor_R35_port_h010_netgen_ascii.msh" \
  "$PROJECT_ROOT/runs/em-10w/Q_Pinc10W_cells.csv" \
  --power 10 --tend 300 --dt 1 \
  --omega106 0.0036 --omega107 0.0036 \
  --damage-threshold 0.98 --save-every 15
```

The unified solver supports power scaling as `Q(P) = (P/10) Q(10 W)` without resolving the electromagnetic model. This assumes linear electromagnetics and temperature-independent properties. Thermal boundaries are adiabatic and no explicit blood vessels are modeled.

## Archive scope and reproducibility

All six supplied research files are preserved byte-for-byte. The README and checksum manifest were added for repository navigation and provenance. No numerical simulations were rerun as part of publication.

The PDF describes a larger supplementary archive, including build files, separate earlier thermal targets, plotting/export scripts, sweep drivers, and result tables. Those additional files are not included here. The thermal source supplied with this repository is `bioheat_R35_unified.cpp`; follow its interface rather than assuming that the older target names in the PDF exist in this upload.

To verify the archived research files on macOS or Linux:

```bash
shasum -a 256 -c SHA256SUMS
```

## Attribution and reuse

Please attribute this work to Viktoria Đorđević and cite the October 2026 internship manuscript, *Three-Dimensional Finite-Element Modeling of Microwave Ablation of a Patient-Derived Liver Tumor Using a Dual-Slot Coaxial Antenna*, Institute of Physics Belgrade.

No software or data reuse license has been selected for this repository. Dataset provenance and references are documented in the manuscript; publication here does not assign new rights to third-party source material.
