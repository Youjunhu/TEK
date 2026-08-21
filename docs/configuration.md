# TEK input and configuration reference

TEK reads a Fortran namelist file named **`input.nmlt`** from the current working directory. The file is read several times by different modules, but all reads use the same file. Every path in the namelists is therefore relative to the directory from which `TEK` is launched.

The checked-in root `input.nmlt` is an ITPA-EP-style electromagnetic three-species example. `itpa_ep/input0.nmlt` is a smaller two-species case, and `itpa_ep/input_tae.nmlt` is a template in which `tfxxx.txt` is replaced by a selected fast-ion temperature profile by the cluster scripts.

## 1. `&magnetic_configuration`

| Parameter | Type / unit | Meaning |
|---|---|---|
| `gfile_name` | character | G-EQDSK-like equilibrium file. The parser expects the standard header, `nx`/`nz` dimensions, flux/profile arrays, safety factor, LCFS, and limiter sections. |
| `reverse_tf` | logical | Reverse the toroidal-field sign convention when required by the input equilibrium. |
| `reverse_ip` | logical | Reverse the plasma-current/poloidal-field sign convention when required by the input equilibrium. |

The equilibrium reader is in `equilibrium.f90`. It uses the file’s `nx` as the one-dimensional equilibrium/radial profile count and reads the two-dimensional poloidal flux array with dimensions `nx` by `nz`.

## 2. `&control_nmlt`

| Parameter | Type / unit | Meaning and constraints |
|---|---|---|
| `pfn_inner` | real, normalized flux | Inner normalized poloidal-flux boundary of the simulation region. |
| `pfn_bdry` | real, normalized flux | Outer normalized poloidal-flux boundary. Must be greater than `pfn_inner`. |
| `nrad` | integer | Number of radial simulation grid points/surfaces. The Poisson and transform code use interior size `nrad - 2`, so use a value comfortably above 2. |
| `mtor` | integer | Number of toroidal grid points in the simulated toroidal segment. Must be even for the FFT conventions used by the code. |
| `mpol` | integer | Number of equilibrium poloidal points. Expected to be odd; `zgrid` spans `[-pi, pi]` with the endpoint duplicated. |
| `ntube` | integer | Number of tube groups in the MPI decomposition. Must divide `numprocs`. |
| `poloidal_angle_type` | character | Magnetic-surface angle construction. Source cases use `straight-field-line`, `equal-volume`, and comments also identify `equal-arc` and `Boozer`. Verify the selected path in `mapping.f90`/`magnetic_coordinates.f90` for a new case. |
| `nsegment` | integer | Toroidal periodicity reduction. The simulated toroidal range is `2*pi/nsegment`; the full torus is represented by the selected segment/harmonics. |
| `nh_min`, `nh_max` | non-negative integers | Toroidal harmonic range used by the field solvers/diagnostics. The control module documents the included harmonics as positive and negative `nsegment` multiples over this range. |
| `dt_omega_i_axis` | real, normalized | Time-step control expressed relative to the reference ion cyclotron angular frequency at the magnetic axis. `main.f90` converts it to seconds and to species-specific normalized steps. |
| `kstart` | integer | Starting iteration. `0` loads fresh markers/fields; a positive value reads `myidxxxxx.pd` restart files and must equal the iteration stored in those files. |
| `kend` | integer | Final iteration. The main loop runs from `kstart + 1` through `kend`. |
| `filter_radial` | logical | Enables the radial filtering path in the Poisson solve. |
| `ismooth` | integer | Selects the smoothing option used by the field solver/filter modules. `0` is the common sample value. Inspect `filter.f90` if using another value. |
| `space_charge_switch` | integer | Switch passed to polarization/field logic for space-charge-related options. The source contains the hook, but the active effect is case-dependent; validate it against the current implementation before relying on it for a study. |
| `diagnosis` | logical | Enables extra geometry/grid diagnostic output and startup visualization routines. It is not the same as the always-open evolution streams. |
| `iplot_mode_structure` | integer | Iteration interval for mode-structure snapshots. The main loop starts these snapshots only after iteration 20,000 and only on the designated tube rank. |
| `store_restart_data` | logical | Writes per-rank unformatted restart files after the run. |
| `fk_switch` | integer | `0` for the GK-only path; `1` initializes, loads, advances, and deposits fully kinetic ions. |
| `adiabatic_electrons` | logical | `true` adds the adiabatic-electron response to the Poisson system and skips the Ampère branch; `false` uses the electromagnetic `A_parallel` path. |
| `nsm` | integer | Number of gyrokinetic species. Every array-valued entry in `&gk_nmlt` must provide exactly `nsm` values unless the intended Fortran namelist behavior is explicitly understood. |

### MPI/grid consistency

For `numprocs` MPI ranks, the executable enforces:

```text
numprocs % ntube == 0
(mpol - 1) % (numprocs / ntube) == 0
```

The root sample has `ntube = 16` and `mpol = 129`, so a valid process count must be divisible by 16 and must make `128 / (numprocs/16)` an integer. The root case’s documented launch of 512 ranks satisfies this (`mpol2 = 32`).

## 3. `&gk_nmlt`

All arrays below are dimensioned by `nsm`.

| Parameter | Type / unit | Meaning |
|---|---|---|
| `mass_gk(:)` | real, kg | Species masses. |
| `charge_gk(:)` | real, elementary-charge units | Charge sign/magnitude before the code multiplies by `elementary_charge`; use `1` for singly charged ions and `-1` for electrons. |
| `nm_gk_per_cell(:)` | integer | Initial markers per `(radial, perturbed-poloidal, toroidal)` cell. Total species markers are `nm_gk_per_cell * nrad * mpol2 * mtor`. |
| `gk_flr(:)` | logical | Whether the species uses finite-Larmor-radius gyro-ring treatment. A zero-FLR electron is represented with `.false.`. |
| `density_file(:)` | character | Two-column radial profile file: coordinate followed by density value. |
| `density_unit(:)` | real | Multiplier applied to the file values before conversion to the code normalization. |
| `density_radcor(:)` | character | Coordinate convention in the profile file. Implemented forms include `toroidal-flux-sqrt`, `poloidal-flux`, and `poloidal-flux-sqrt`. |
| `temperature_file(:)` | character | Two-column temperature profile file. |
| `temperature_unit(:)` | real | Multiplier applied to file temperatures. The sample uses `1.6022e-16` to convert keV-like values to joules where appropriate. |
| `temperature_radcor(:)` | character | Radial-coordinate convention for the temperature file. |
| `gk_spatial_loading_scheme(:)` | integer | `1`: uniform loading in magnetic `(psi, theta, alpha)` coordinates; `2`: uniform real-space loading. |
| `gk_velocity_loading_scheme(:)` | integer | Velocity initialization selector. Sample cases use `2`, the isotropic/Maxwellian-oriented path. |
| `gk_nonlinear(:)` | integer | Per-species switch multiplying nonlinear field terms in `push_gk_weight.f90`. Use one value per species. |

The profile loader reads at most 3000 rows and requires more than one row. Profile values are linearly interpolated onto the simulation’s radial grid. Files must contain two whitespace-separated columns and must cover the requested coordinate range adequately; extrapolation behavior is not a substitute for a physically valid profile.

The code establishes the shared normalization from species 1. In particular, `qu`, `tu`, `vu`, and `nu` are assigned from the first species, and the main time step is based on its `dtao_gk(1)`. Choose the first species intentionally.

## 4. `&adiabatic_electron_nmlt`

This group is read only when `adiabatic_electrons = .true.`.

| Parameter | Type / unit | Meaning |
|---|---|---|
| `density_file` | character | Two-column electron density profile. |
| `density_unit` | real | Density scale multiplier; the source comments expect SI `m^-3` values after scaling. |
| `density_radcor` | character | Profile coordinate convention. |
| `temperature_file` | character | Two-column electron temperature profile. |
| `temperature_unit` | real | Temperature scale multiplier; the loader converts into the code’s energy normalization. |
| `temperature_radcor` | character | Temperature profile coordinate convention. |

The resulting adiabatic response is added to the Poisson matrix diagonal for the zero toroidal harmonic. The electromagnetic Ampère solve is skipped in this mode.

## 5. `&fk_nmlt`

This group is read by `initialize_fk()` only when `fk_switch = 1`.

| Parameter | Type / unit | Meaning |
|---|---|---|
| `mass_i` | real, kg | Fully kinetic ion mass. |
| `charge_i` | real, C | Fully kinetic ion charge in coulombs (unlike `charge_gk`, which is specified in elementary-charge units). |
| `ti0` | real, keV | Reference FK ion temperature used for velocity/loading normalization. |
| `ni0` | real, `m^-3` | Reference FK ion density. |
| `kappa_ni` | real, `Ln^-1` | Density-gradient parameter used by the FK weight model. |
| `kappa_ti` | real, `Ln^-1` | Temperature-gradient parameter used by the FK weight model. |
| `nmarker_i_per_cell` | integer | Total initial FK ion markers per simulation cell. The total is `nmarker_i_per_cell * nrad * mpol2 * mtor`, divided across MPI ranks. |
| `ion_spatial_loading_scheme` | integer | `1`: magnetic-coordinate loading; `2`: real-space loading. |
| `ion_velocity_loading_scheme` | integer | `1`: Cartesian velocity-coordinate loading; `2`: isotropic Gaussian loading. |
| `fk_nonlinear` | integer | Switch used by the FK weight/orbit path for nonlinear terms. |

The FK arrays are allocated at a fixed over-capacity derived from the initial marker count. FK markers are advanced in cylindrical coordinates and repeatedly converted to magnetic coordinates for sorting/deposition.

## 6. Profile file format

A profile is plain text with two columns and no header. For example:

```text
0.0000000000000000  5.8288292471993057E+019
0.0078125000000000  5.7815345306430030E+019
```

The first column is interpreted according to the selected `*_radcor` string. The second column is multiplied by `*_unit`. The checked-in `cbc/*.txt` files are poloidal-flux examples; `itpa_ep/*.txt` files are mainly toroidal-flux-sqrt examples.

## 7. Recommended configuration procedure

1. Copy a nearby case file rather than editing the root example in place.
2. Set the equilibrium and verify the G-EQDSK path from the intended working directory.
3. Set `nsm` and provide every `nsm`-length GK value explicitly.
4. Choose `numprocs` and `ntube` together, then check both decomposition constraints.
5. Confirm that `mtor` is even and `mpol` is odd.
6. Check the profile coordinate convention and units for every species.
7. Start with a short run (`kend` only a few time steps beyond `kstart`) and inspect startup diagnostics, `xgrid.txt`, profile outputs, and field evolution before launching a production run.
8. Preserve the exact input file and MPI layout with any restart files.
