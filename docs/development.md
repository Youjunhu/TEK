# Source tree and development guide

## 1. Source inventory

The `makefile` source order is the authoritative compilation order. The table below groups the source files by responsibility.

### Foundations and numerical utilities

| File | Main units / role |
|---|---|
| `modules.f90` | `constants`, `normalizing`, equilibrium/grid state modules, MPI decomposition, controls, and perturbation-field storage. |
| `math.f90` | Interpolation module plus geometry, derivatives, complex solves, Fourier-mode solver wrappers, and numerical helpers. |
| `splines.f90` | LAPACK declarations/wrappers, text utilities, cubic/2-D spline routines, indexing and polynomial helpers. |
| `pputil_yj.f90` | MPI particle movement/sorting utilities used to redistribute marker arrays. |
| `transform.f90` | FFTW3 plans, toroidal/radial Fourier transforms, sine transforms, and retained Numerical Recipes transforms. |

### Equilibrium and coordinates

| File | Main units / role |
|---|---|
| `equilibrium.f90` | G-EQDSK reader, poloidal-flux interpolation, magnetic-field functions, equilibrium-derived profiles. |
| `magnetic_coordinates.f90` | Magnetic surfaces, angle construction, metrics, Jacobians, toroidal shift, drift tables, and coordinate helpers. |
| `mapping.f90` | Cylindrical/magnetic mapping table construction and interpolation. |
| `profiles.f90` | Generic radial profile type plus adiabatic-electron and GK profile objects/functions. |

### GK particles and fields

| File | Main units / role |
|---|---|
| `gk_module.f90` | Species arrays, GK initialization, marker capacity, profile loading, and GK marker sorting. |
| `load_gk.f90` | Marker loading, equilibrium distribution functions, velocity sampling, and particle-to-guiding-center helpers. |
| `gyro_ring.f90` | Four-point gyro-ring construction and gyro phase initialization. |
| `gyro_average.f90` | Interpolation/averaging of perturbation fields at gyro points. |
| `drift.f90` | Guiding-center drift and mirror-force evaluation from geometry tables and fields. |
| `push_gk_orbit.f90` | Guiding-center characteristic push and periodic coordinate wrapping. |
| `push_gk_weight.f90` | GK delta-f weight equation for Maxwellian/slowing-down variants and nonlinear switches. |
| `deposit_gk.f90` | GK density and parallel-current deposition with trilinear/inter-cell weighting. |
| `gk_polarization.f90` | Polarization matrix construction, including the active Padé approximation and experimental direct forms. |
| `poisson.f90` | Matrix assembly/factorization and harmonic-by-harmonic Poisson solve. |
| `ampere.f90` | Electromagnetic parallel-vector-potential solve, split-field evolution, derivatives, and pullback. |

### Optional FK path

| File | Main units / role |
|---|---|
| `fk_module.f90` | FK ion arrays and initialization namelist. |
| `boris.f90` | Cylindrical full-orbit Boris pusher, normalization, and GK/FK conversion helpers. |
| `force.f90` | Perturbed fields and electron/full-orbit force functions. |
| `push_fk_orbit.f90` | FK cylindrical-to-magnetic coordinate transform, loss handling, and cleanup. |
| `push_fk_weight.f90` | FK density/temperature functions, ion weight push, and first-step wrapper. |
| `deposit_fk.f90` | FK density/current deposition. |

### Communication, filtering, diagnostics, and persistence

| File | Main units / role |
|---|---|
| `communication_connection.f90` | MPI boundary/cut communication and source merging. |
| `derivatives_in_xyz.f90` | Distributed x/y/field-line derivatives. |
| `filter.f90` | Toroidal/radial filtering, sine filtering, reconstruction, and smoothing. |
| `diagnosis.f90` | Evolution streams, spectra, mode structures, flux/entropy checks, grid and field-line diagnostics. |
| `restart.f90` | Per-rank unformatted checkpoint writer/reader. |
| `main.f90` | Executable orchestration and `&control_nmlt` reader. |

## 2. Build ordering

The current source order is:

```text
modules
pputil_yj
splines
math
equilibrium
magnetic_coordinates
mapping
transform
fk_module
profiles
gk_module
load_gk
deposit_gk
deposit_fk
communication_connection
filter
derivatives_in_xyz
gyro_ring
gyro_average
gk_polarization
poisson
ampere
force
boris
push_fk_orbit
drift
push_gk_orbit
push_gk_weight
push_fk_weight
diagnosis
restart
main
```

When adding a module, place it after all modules whose `.mod` files it uses. If the separate-object build reports missing module files while the single-file build works, update the source order and/or add explicit make dependencies.

## 3. Shared state conventions

The code uses extensive module-level `save` state rather than passing a complete simulation object. Before changing an array, identify:

- its allocation site;
- its indexing convention and boundary/ghost points;
- its owning communicator/rank;
- whether it is normalized or SI;
- whether it is current-time, midpoint, or old-time state;
- whether `(:,:)` means local ownership or global logical coordinates.

Common examples:

- `potential`, `apara`, and derivatives are `(mtor, nrad, 2)` local left/right field-cell arrays;
- `phi_dft` and `apara_dft` are `(0:mtor-1, nrad-2)` complex harmonic arrays;
- GK marker arrays are `(nmmax, nsm)` with actual active counts in `nm_gk(:)`;
- `x_ring`, `y_ring`, `z_ring` are `(gyro_npt, nmmax, nsm)`;
- `xgrid` is radial, `ygrid` toroidal, and `zgrid` poloidal despite the source’s historical x/y/z naming.

## 4. Units and normalization

`constants` defines SI constants and `p_ = kind(1.0d0)`. `normalizing` fixes `Ln = 1 m` and `bn = 1 T` and stores the runtime species-1 normalizations. Temperatures are commonly represented as keV in profile routines and converted with `kev`; charges in GK input are elementary-charge units and are converted to coulombs, while FK `charge_i` is already in coulombs.

When modifying equations, trace every factor involving:

- `qu`, `tu`, `vu`, `nu`;
- `dtao_gk`, `dtao_main`, `dtao_fk`;
- `w_unit` and `ps_vol_*`;
- `kev`, `elementary_charge`, `epsilon0`, `mu0`, and `c`.

The first GK species defines shared normalization, so changing species ordering can change numerical scaling even if the physical species set is unchanged.

## 5. Debugging workflow

1. Build the single-file executable with bounds checking and backtraces.
2. Run a short case with one or a few valid MPI layouts.
3. Compare startup diagnostics (`q1.txt`, `xgrid.txt`, `profiles_ns*.txt`, geometry files) against the equilibrium/profile input.
4. Check marker counts and `nmmax` headroom after sorting.
5. Inspect the first `phi_evolution*`, `apara_evolution*`, and `bperp_evolution*` files for NaNs, zeros, or unexpected growth.
6. Use the existing geometry/field-line diagnostics in a scratch directory.
7. Only after a serializable short run passes should you test restart or change MPI decomposition.

The default compiler flags in the fallback branch already include `-fbounds-check`, `-fimplicit-none`, and several warnings. Add `-fbacktrace`, `-ffpe-trap=invalid,zero,overflow`, or optimization controls locally as appropriate for the compiler and case.

## 6. Known limitations and maintenance risks

These are observations from the current source, not a claim that the scientific model is invalid:

- There is no automated test suite or CI configuration in the repository.
- The makefile contains hard-coded paths and historical cluster commands; the Linux fallback’s absolute library filenames may need overrides.
- The code writes many fixed diagnostic filenames and can overwrite files in the working directory.
- Restart files are compiler-/layout-dependent Fortran unformatted records and save GK state/fields, while the FK state is not included in the active record.
- Most configuration arrays are allocated but not default-initialized before namelist reading. Every `nsm`-length value should therefore be written explicitly. For example, a scalar `gk_nonlinear = 1` in a two-species case leaves the second array element unspecified under normal Fortran namelist semantics.
- Some control variables and diagnostic routines are legacy or experimental. A parameter appearing in a namelist is not by itself proof that every possible value has a complete, tested implementation.
- The source contains old-style LAPACK declarations and floating-point equality tests; these produce warnings with modern compilers.
- The default time loop writes diagnostics at fixed iteration intervals and only after a hard-coded mode-structure threshold (`kt > 20000`).
- There are no schema/version checks for profile files, equilibrium files, output files, or restart records.

## 7. Adding a new case

Keep case assets separate from the executable source when possible:

```text
my_case/
  input.nmlt
  equilibrium.gfile
  profiles/
    ni.txt
    ne.txt
    ti.txt
    te.txt
  run.sh
  README.md
```

Use paths relative to the launch directory, or launch from the case directory after copying/symlinking the executable. Record the compiler command, Git revision, MPI process count, environment modules, and all output filenames. For a production benchmark, retain both the original profile files and the generated/modified input rather than relying on a shell substitution script to reconstruct them.
