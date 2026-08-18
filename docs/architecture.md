# TEK architecture and numerical workflow

## 1. Purpose and model

TEK advances a delta-f representation of low-frequency tokamak turbulence with a particle-in-cell method. The primary path is gyrokinetic (GK): marker guiding centers carry magnetic coordinates, parallel velocity, magnetic moment, and a perturbation weight. The code supports multiple GK species through `nsm`, including a common configuration with finite-Larmor-radius ions and zero-FLR electrons. When `adiabatic_electrons = .false.`, the electromagnetic parallel vector potential is evolved by the mixed-variable/Ampère path; when it is true, the electron response is supplied through an adiabatic-electron density/temperature model in the Poisson matrix.

A separate optional path (`fk_switch = 1`) initializes and advances fully kinetic ions with cylindrical full-orbit variables and a Boris pusher. The FK path is coupled to the field deposition/weight machinery, but the restart writer currently serializes GK state and fields only; see the limitations section.

The implementation is organized around a structured magnetic-coordinate grid:

- radial coordinate `x`/`radcor`, generated from `pfn_inner` to `pfn_bdry`;
- poloidal coordinate `z`/`theta`, with `mpol` equilibrium points and a distributed perturbed-field grid of `mpol2 = numprocs / ntube` points per process group;
- toroidal coordinate `y`/`alpha`/toroidal angle, with a computational toroidal range `2*pi/nsegment` sampled by `mtor` points.

The equilibrium itself is read in cylindrical `(R,Z)` coordinates from a G-EQDSK-like file and transformed into magnetic coordinates before marker loading.

## 2. Startup sequence

`main.f90` is the executable entry point. Its startup sequence is:

1. `MPI_INIT`, followed by `MPI_COMM_SIZE` and `MPI_COMM_RANK`.
2. Read `input.nmlt` through `read_parameters()` (`&control_nmlt`).
3. Derive the decomposition:
   - `nvp = numprocs / ntube`;
   - `GCLR = myid / ntube` identifies the poloidal grid rank;
   - `TCLR = mod(myid, ntube)` identifies the tube/toroidal rank;
   - `GRID_COMM` and `TUBE_COMM` are created with `MPI_COMM_SPLIT`.
4. Validate `numprocs % ntube == 0` and `(mpol - 1) % (numprocs / ntube) == 0`.
5. Read/process the equilibrium, construct magnetic surfaces and metrics, build the cylindrical-to-magnetic mapping table, and prepare geometry tables for drift calculations.
6. Initialize GK profiles and marker arrays. The first species establishes the shared normalization (`qu`, `tu`, `vu`, `nu`, and `dtao_main`).
7. Allocate perturbation fields. If `kstart == 0`, fields are zeroed and markers are loaded from equilibrium distributions; if `kstart > 0`, a per-rank binary restart is read.
8. Set gyro phases and, if enabled, initialize/load the FK ion population.
9. Initialize FFTW plans, open evolution output files, and factorize the Poisson and Ampère matrices.
10. Enter the time loop.

The code writes a substantial amount of geometry/profile information during startup. This is useful for validation but means that a successful process launch is not equivalent to a physically validated case.

## 3. Geometry and coordinate construction

### 3.1 Equilibrium input

`equilibrium.f90` reads the magnetic configuration namelist and then parses the selected G-EQDSK file. It loads the poloidal flux mesh, axis/LCFS fluxes, pressure/current profiles, safety factor, LCFS, and limiter data. The poloidal flux is interpolated with the interpolation/spline utilities, and normalized flux tables are generated.

`reverse_tf` and `reverse_ip` are passed from the namelist into the equilibrium reader to control sign conventions for toroidal field/current orientation. Use them only when the equilibrium convention requires it; changing either changes the signs used by downstream geometry and field routines.

### 3.2 Magnetic surfaces

`magnetic_coordinates.f90` selects `nrad` surfaces between `pfn_inner` and `pfn_bdry`, finds contours, constructs the requested poloidal angle (`straight-field-line`, `equal-arc`, `equal-volume`, or `Boozer` where implemented), and computes:

- `r_mc`, `z_mc` surface coordinates;
- grid arrays `xgrid`, `ygrid`, `zgrid`;
- Jacobian and metric tensors;
- gradients and dot products of `psi`, `theta`, and generalized toroidal angle `alpha`;
- magnetic field and curvature/drift lookup tables;
- computational volume and radial cell volumes.

The radial grid is uniform in the selected normalized flux interval before conversion to the internal radial coordinate. The code writes `xgrid.txt` and several geometry diagnostics on rank zero.

### 3.3 Mapping

`mapping.f90` builds a lookup table for converting cylindrical marker positions to magnetic coordinates. This is used by the fully kinetic path, which advances in cylindrical coordinates but must sort and deposit markers in magnetic coordinates. The reverse interpolation path in `map_to_mc` is used by marker loading and geometry operations.

## 4. MPI decomposition and communication

The implementation uses a two-dimensional logical organization of MPI ranks:

- `GCLR` (grid communicator coordinate) partitions the poloidal/field-line direction;
- `TCLR` (tube communicator coordinate) groups ranks that share the same poloidal location across the toroidal decomposition.

The solver assumes periodicity in toroidal direction and performs explicit neighbor exchange across the poloidal cut. `communication_connection.f90` handles:

- copying/connecting scalar values across the theta cut;
- merging left/right source contributions from marker deposition;
- updating field derivatives at a process boundary;
- field-value exchange between neighboring cells;
- retrieval of nearby values along field lines.

The decomposition constraints are structural:

```text
numprocs % ntube == 0
(mpol - 1) % (numprocs / ntube) == 0
```

The `mpol2 = numprocs/ntube` perturbed-field poloidal count is used in array dimensions and in the `dtheta2` spacing. A case that satisfies only the first condition can still fail during startup if the second condition is not met.

## 5. Marker state and loading

### 5.1 GK state

`gk_module.f90` owns per-species arrays. Important state includes:

- `xgc`, `zgc`, `ygc`: guiding-center magnetic coordinates;
- `vpar_gk`: parallel velocity in normalized units;
- `mu_gk`: magnetic moment;
- `v_gk`: speed used by the distribution/weight equations;
- `w_gk`: delta-f marker weights;
- `ps_vol_gk`: phase-space volume/marker normalization;
- `lost_gc`: radial loss flag;
- `_mid` arrays for the intermediate Runge–Kutta stage;
- `x_ring`, `y_ring`, `z_ring`: four-point gyro rings (`gyro_npt = 4`).

The initial total markers for species `s` are:

```text
total_nm_gk(s) = nm_gk_per_cell(s) * nrad * mpol2 * mtor
```

Each rank initially receives `total_nm_gk(s)/numprocs`. Arrays are over-allocated with a factor of 2.5 relative to the maximum initial per-rank population to allow for redistribution.

`load_gk.f90` supports the configured spatial and velocity loading schemes, initializes Maxwellian-like distributions, assigns marker weights, and provides distribution/profile functions used later by the weight pusher. `sort_gk_markers` redistributes markers according to poloidal location.

### 5.2 Gyro ring and field sampling

`gyro_ring.f90` computes four gyro-phase points for each guiding center. `gyro_average.f90` interpolates perturbation fields to those ring points and averages them back to the guiding center. The averaged fields feed `drift.f90`, which calculates radial/poloidal/toroidal drifts and mirror-force terms.

### 5.3 FK state

When `fk_switch = 1`, `fk_module.f90` allocates cylindrical position/velocity arrays and magnetic-coordinate work arrays. `boris.f90` implements cylindrical full-orbit Boris operations, coordinate normalization, and particle/guiding-center conversion helpers. `push_fk_orbit.f90` maps FK markers to magnetic coordinates and removes/out-of-domain markers. `push_fk_weight.f90` evolves ion delta-f weights, while `deposit_fk.f90` deposits FK density/current contributions.

## 6. One time step

The main loop is a two-stage second-order Runge–Kutta update. In simplified form, each iteration does the following.

### Stage 1: half step

1. If the FK path is active, advance the FK particles by the first half-step and deposit their source.
2. For each GK species:
   - gyro-average fields at current marker positions;
   - compute guiding-center drifts and mirror force;
   - advance weights by half a GK time step;
   - advance guiding centers by half a GK time step;
   - mark radial losses and redistribute markers;
   - build gyro rings and deposit density and parallel current.
3. Solve Poisson for the half-step electrostatic potential.
4. For electromagnetic runs, evolve the split scalar part of `A_parallel` and solve Ampère for the half-step vector potential.

### Stage 2: full step

1. Advance the FK path through its second stage and deposit its source, if enabled.
2. For each GK species, gyro-average the intermediate fields, recompute drifts, advance weights and guiding centers through a full step, sort/mark particles, form gyro rings, and deposit again.
3. Solve Poisson again.
4. In electromagnetic mode, update `A_parallel`, solve Ampère, and perform the mixed-variable weight pullback/resplit.
5. Write evolution diagnostics at their configured cadence and optionally write mode-structure snapshots.

The relevant high-level calls are visible in `main.f90` around the loop beginning after matrix preparation. The actual equations and normalization factors are implemented in the pusher, drift, deposition, polarization, Poisson, and Ampère modules.

## 7. Field solves

### 7.1 Poisson

`gk_polarization.f90` constructs a species polarization matrix. The active implementation uses a finite-Larmor-radius response with a Padé-style approximation for the gyroaveraging factor; alternative direct numerical polarization builders remain in the source as experimental routines. `poisson.f90` sums species matrices, optionally adds the adiabatic-electron diagonal response for the `n=0` component, and factorizes each toroidal harmonic with LAPACK `ZGETRF`.

At runtime, marker density is merged across the left/right process boundaries, normalized by cell volume and the reference density, transformed into toroidal Fourier space, solved harmonic-by-harmonic, and transformed back. Radial filtering/smoothing and zero-mode handling are available through `filter.f90`, controlled by `filter_radial` and `ismooth`.

The allocated field arrays include both sides of a local poloidal cell (`(:,:,1)` and `(:,:,2)`) so that derivatives and communication can be performed across neighboring field-line cells.

### 7.2 Ampère / electromagnetic response

When `adiabatic_electrons = .false.`, `ampere.f90` prepares and solves the parallel-vector-potential system. The implementation separates `A_parallel` into `A_s` and `A_h` for the mixed-variable pullback method:

- `apara_s` is evolved from the electrostatic parallel derivative;
- `apara_h` is obtained from the current/Ampère solve;
- `apara` is the reconstructed total field;
- the weight pullback updates GK weights after the full step.

The module also computes derivatives, a Laplacian, and skin-current residual diagnostics. In adiabatic-electron mode, the Ampère branch in `main.f90` is skipped.

### 7.3 Spectral transforms

`transform.f90` wraps FFTW 3 plans for toroidal/radial transforms and a radial sine transform, plus retained Numerical Recipes-style routines for compatibility/testing. The toroidal grid is periodic; the radial sine transform operates on the interior `nrad - 2` points. FFTW plans are initialized once after the geometry and field dimensions are known.

## 8. Diagnostics and restart

`diagnosis.f90` contains field evolution, harmonic spectra, mode-structure output, particle/heat flux, entropy, grid/field-line checks, and geometry visualization routines. `main.f90` always opens the primary evolution streams; `diagnosis` enables additional geometry/testing output during startup, while `iplot_mode_structure` controls periodic poloidal-plane snapshots after the configured iteration threshold.

`restart.f90` writes one unformatted binary file per MPI rank:

```text
myid00000.pd, myid00001.pd, ...
```

The file contains the saved iteration, GK marker counts and arrays, loss flags, and electrostatic/electromagnetic fields. Restarting requires the same `kstart` value as the saved iteration and the same decomposition/array layout. The reader reconstructs `apara_s` from `apara` and resets the complementary split fields. FK particle arrays are currently commented out of the restart record.
