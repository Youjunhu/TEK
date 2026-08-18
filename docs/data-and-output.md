# TEK data files, runtime outputs, and diagnostics

## 1. Equilibrium files

The magnetic reader expects a G-EQDSK-like text file. The checked-in examples include:

- `cbc/gfile_circular` — circular benchmark equilibrium;
- `itpa_ep/gfile_itpa_ep` — high-resolution ITPA-EP equilibrium;
- `itpa_ep/gfile_itpa_ep129x129`, `gfile_itpa_ep257x257`, `gfile_itpa_ep513x513` — different equilibrium resolutions;
- `itpa_ep/gfile_itpa_ep_GTAW_output` — another saved equilibrium variant.

The file contains a header with `nx`, `nz`, geometry scalars, axis/LCFS fluxes, toroidal-field/current profiles, the 2-D poloidal flux, safety factor, LCFS coordinates, and limiter coordinates. `equilibrium.f90` reads these sections with fixed formats. A file that is merely “similar to” an EQDSK but has a different record layout is not necessarily compatible.

## 2. Radial profiles

All profile files are two-column whitespace-delimited text without a header:

```text
radial_coordinate  value
```

The loader supports these coordinate labels:

- `toroidal-flux-sqrt` — first column is square-root toroidal flux coordinate;
- `poloidal-flux` — first column is normalized poloidal flux;
- `poloidal-flux-sqrt` — first column is square-root normalized poloidal flux.

`cbc/` contains Cyclone-Base-Case-style density/temperature examples. `itpa_ep/` contains density/temperature profiles for thermal and fast-ion examples. `itpa_ep/profiles.py` generates several of the flat or parameterized files and `itpa_ep/tae.py` prints simple TAE scale estimates.

The profile reader accepts up to 3000 rows, requires at least two, scales values by the configured unit, and linearly interpolates onto the simulation radial grid. Keep profile files in the working directory tree addressed by `input.nmlt`.

## 3. Startup/geometry output

The following files are written by the normal initialization path, usually by rank zero unless noted otherwise. Existing files may be overwritten.

| File | Producer / content |
|---|---|
| `q1.txt` | Equilibrium normalized flux, toroidal-flux coordinate, and safety-factor table. |
| `xgrid.txt` | Simulation radial coordinate and toroidal-flux coordinate. |
| `profiles_ns1.txt`, `profiles_ns2.txt`, ... | Interpolated species profiles and radial derivatives. |
| `mapping_table.txt` | Cylindrical-to-magnetic mapping data. |
| `rzgrid.txt` | Cylindrical grid used for mapping diagnostics. |
| `grid.txt`, `grid_z*.txt`, `inner_bdry.txt` | Grid/field-line visualization diagnostics from `diagnosis.f90`. |
| `mag_surf_shape*.txt`, `minor_r.txt`, `theta_line.txt` | Magnetic-surface and poloidal-coordinate diagnostics. |
| `grad_theta_alpha.txt`, `qhat.txt`, `gradxz.txt`, `grady.txt` | Metric/gradient diagnostics. |
| `rectangular`, `limiter.txt`, `lcfs.txt`, `lcfs2.txt` | Boundary visualization/export files produced by geometry routines. |

The exact set depends on `diagnosis`, selected angle type, and which diagnostic routines are enabled. Some diagnostics are available as callable routines but are commented out in the main program.

## 4. Main evolution streams

`main.f90` opens files using a nine-digit, zero-padded `kstart` suffix. The designated evolution-writing rank is `numprocs/2`, corresponding to the code’s theta-zero location; rank zero writes the `B_perp` stream.

| File pattern | Content / cadence |
|---|---|
| `phi_evolution000000000.txt` | Real-space potential mode/evolution diagnostic. Written every 10 iterations by the designated rank. |
| `apara_evolution000000000.txt` | Real-space parallel-vector-potential evolution diagnostic. Written every 10 iterations. |
| `phi_n_evolution000000000.txt` | Toroidal harmonic data for `phi_dft`. Written every 10 iterations. |
| `apara_n_evolution000000000.txt` | Toroidal harmonic data for `apara_dft`. Written every 10 iterations. |
| `bperp_evolution000000000.txt` | Perpendicular magnetic perturbation diagnostic, written by rank zero every 100 iterations. |

The suffix is the starting iteration, not the final iteration. A restart therefore creates a new stream suffix for the restarted `kstart`.

`mode_evolution` and `nharmonic_evolution` live in `diagnosis.f90`; their exact column layout is source-defined rather than documented by a schema file. Preserve the corresponding input, source revision, and grid sizes with any analysis product.

## 5. Mode-structure output

When `TCLR == 0`, `kt > 20000`, and `mod(kt-1, iplot_mode_structure) == 0`, the main loop writes poloidal-plane mode-structure files through:

```fortran
call mode_structure_in_poloidal_plane(kt, potential, 'Phi')
call mode_structure_in_poloidal_plane(kt, apara,     'Apara')
```

Other XY/XZ/YZ and harmonic output routines are present in `diagnosis.f90` but are commented out in the default loop. The resulting files are intended for visualization and are not a stable public file format.

## 6. Restart files

With `store_restart_data = .true.`, each rank writes an unformatted Fortran binary file:

```text
myid00000.pd
myid00001.pd
...
```

The record contains, in order:

1. saved iteration `kend`;
2. `nm_gk`;
3. GK weights, phase-space volume, loss flags, speed, coordinates, magnetic moment, and parallel velocity;
4. `potential`, `phix`, `phiy`, `phiz`;
5. `apara`, `ax`, `ay`, `az`.

Because this is an unformatted compiler-dependent record with assumed-shape/allocated-array state, it should be treated as a same-build/same-layout checkpoint, not a portable interchange format. The restart reader requires exact iteration equality and reconstructs part of the mixed-variable state after reading.

## 7. Diagnostic routines

`diagnosis.f90` includes more capabilities than the default main loop invokes:

- particle and heat flux;
- entropy and lost-marker counts;
- parallel current and weight checks;
- grid and surface drawings;
- toroidal harmonic/mode evolution;
- perpendicular magnetic perturbation;
- XY/XZ/YZ/poloidal mode structures;
- field-line tracing and Poincaré-like outputs;
- safety-factor and boundary-touch checks.

These routines often write fixed filenames and are best run in a scratch copy of a case. They may assume rank ownership, communicator state, or arrays initialized by the normal startup sequence.

## 8. Output hygiene

A production directory can contain large text diagnostics and per-rank restart files. Recommended practice:

- run each case in a dedicated directory;
- copy/record the exact `input.nmlt`, equilibrium, and profile files;
- preserve stdout/stderr and the executable build command;
- do not mix restart files from different MPI layouts;
- move diagnostic files before relaunching if overwrite would destroy a previous result;
- use the `itpa_ep` scripts only on the intended cluster, after reviewing their `rsync`, scheduler, and module commands.
