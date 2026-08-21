# Build, run, restart, and benchmark workflows

## 1. Build prerequisites

The source is Fortran 90/2003-style code and the executable links against:

- MPI Fortran (`mpif90` or site-specific equivalent);
- FFTW 3 (`fftw3.f03` plus the FFTW library);
- BLAS;
- LAPACK.

The code uses FFTW’s C binding from Fortran in `transform.f90` and LAPACK routines such as `ZGETRF`/`ZGETRS` in the field solvers and spline utilities.

## 2. Build system

The legacy `makefile` provides two build styles:

### Separate-object build

The default target `TEK` compiles the ordered source list into `.o` files and links them. The pattern rule has `modules.f90` as a prerequisite for every source, but Fortran module dependency ordering is primarily encoded by the order in `f90source`.

```bash
make
```

### Single-file build

`merge_to_one_file` concatenates all sources in dependency order into `TEK_in_one_file.f90`. `compile_one_file` then compiles and links that file as one translation unit.

```bash
make compile_one_file
```

The single-file target is often the most reproducible option for this repository because it avoids stale `.mod`/`.o` dependency issues and mirrors the order explicitly checked into the makefile.

## 3. Portable Linux build

The default Linux branch in the checked-in makefile refers to absolute library filenames. On systems where those files are installed through the normal linker search path, override the variables:

```bash
make compile_one_file \
  lapack_location=-llapack \
  blas_location=-lblas \
  fftw_lib=-lfftw3
```

A successful build in the current repository environment was verified with this command. It produced only legacy Fortran warnings (old-style character declarations and floating-point equality comparisons), not compile/link errors.

If `fftw3.f03` is not found, add the include directory explicitly:

```bash
make compile_one_file \
  fftw_include=-I/path/to/fftw/include \
  lapack_location=-llapack \
  blas_location=-lblas \
  fftw_lib=-lfftw3
```

The machine branches (`edison`, `cori`, `tianhe`, `3m`, `hfc`, `sm`, and `thex`) encode historical compiler/module/library paths. They should be treated as site templates, not portable defaults.

## 4. Cleaning

```bash
make clean
```

The clean target removes the executable, object files, Fortran module files, `.pd` files under `ms/`, and other generated build artifacts matching its patterns. Review the target before using it in a directory containing valuable runtime data.

## 5. Launching

Run from the project root (or adjust all input paths accordingly):

```bash
export OMP_NUM_THREADS=1
mpirun -n 16 ./TEK
```

The `run` target uses `mpiexec -n 16` and sets `OMP_NUM_THREADS=1`, but it is only a convenience target. The process count must match the input decomposition rules. For the root sample (`ntube=16`, `mpol=129`), a 16-rank launch gives `mpol2=1` and satisfies the divisibility rule, while the sample production launch uses 512 ranks and gives `mpol2=32`.

Before a long launch, check:

```bash
# rank count must be divisible by ntube
# (mpol - 1) must be divisible by (numprocs / ntube)
```

The executable prints the selected namelists on rank zero, grid/volume information, marker counts, normalizations, time step, beta, and CFL-related estimates. Capture stdout/stderr with the scheduler for reproducibility.

## 6. First-run smoke procedure

There is no automated test suite. Use a manual smoke run:

1. copy an existing case to a scratch directory;
2. reduce `kend` to a few steps;
3. choose a process count satisfying the constraints;
4. build with diagnostic compiler flags (`-fbounds-check`, `-fbacktrace`, and, if supported, floating-point traps);
5. verify that the equilibrium/profile files open;
6. inspect `xgrid.txt`, `q1.txt`, `profiles_ns*.txt`, and the first evolution files;
7. only then increase the iteration count and marker/grid resolution.

The source contains several debug/test routines, but they are not wired into a test runner and many are not part of the normal executable path.

## 7. Restart workflow

To write restart files, set:

```fortran
store_restart_data = .true.
```

At the end of the run, each rank writes a file named `myidxxxxx.pd`, where `xxxxx` is the zero-padded MPI rank. To restart:

1. preserve all per-rank `.pd` files;
2. use the same MPI rank count, `ntube`, grid dimensions, species count, and marker allocation parameters;
3. set `kstart` to exactly the iteration stored in the files;
4. set a later `kend`;
5. launch from the directory containing the restart files and input data.

The restart reader explicitly aborts if the input `kstart` does not equal the stored `kend`. It also restores GK fields/markers but not the commented-out FK particle record. Do not assume that changing decomposition or enabling/disabling FK is restart-compatible.

## 8. Checked-in benchmark/case assets

### Root `input.nmlt`

A three-GK-species electromagnetic configuration using the `itpa_ep/gfile_itpa_ep513x513` equilibrium and ion/fast-ion/electron profiles. The sample has `nsm=3`, `ntube=16`, `mtor=16`, `mpol=129`, and `nrad=258`.

### `itpa_ep/input_tae.nmlt`

A two-ion/electron-style ITPA-EP/TAE template. The scripts replace `tfxxx.txt` with a selected `tf100.txt`, `tf200.txt`, etc. before submission.

### `itpa_ep/input0.nmlt`

A smaller two-species example with a coarser radial interval and restart enabled. Its scalar `gk_nonlinear = 1` assignment should be expanded to one value per species if used as a strict reproducibility input; Fortran namelist assignment of a single value to an array does not initialize unspecified elements.

### Cluster scripts

- `itpa_ep/hfc.sh` builds on the HFC platform, creates a Slurm script, copies a case directory, substitutes a fast-ion temperature profile, and submits 512 MPI ranks.
- `itpa_ep/thex.sh` loads Intel/MPI/FFTW/LAPACK modules, builds on TheX, copies the case, substitutes a profile, and submits through `yhbatch`.

Both scripts assume the source tree and scheduler commands are available on the named cluster. They are not local Linux launch scripts.

## 9. Performance and scaling considerations

Memory grows with:

- `nrad * mpol2 * mtor` for field/source arrays;
- `nmmax * nsm` for GK marker arrays;
- the Poisson/Ampère matrices, whose radial solve size is `nrad - 2` for each included harmonic.

Increasing `nsm`, `nm_gk_per_cell`, `nrad`, or `mtor` can multiply memory and setup time. Increasing MPI ranks changes `mpol2`, marker distribution, communicator topology, and restart file count; it is not a transparent performance-only change.
