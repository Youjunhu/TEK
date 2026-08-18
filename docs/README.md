# TEK documentation

This documentation was written from the source tree and checked-in case files. It describes the executable that is built from `main.f90`, the default input contract, the MPI layout, the numerical stages, the generated files, and the current development limitations.

## Reading order

1. [Architecture and numerical workflow](architecture.md) — what the solver computes and how one time step is assembled.
2. [Input and configuration reference](configuration.md) — every namelist group, parameter, unit, and consistency constraint.
3. [Build and run](build-and-run.md) — dependencies, portable compilation, MPI launch, restart procedure, and benchmark cases.
4. [Data and output](data-and-output.md) — input file formats, runtime files, diagnostics, and restart artifacts.
5. [Source and development guide](development.md) — file-by-file source map, build ordering, debugging, and known limitations.

## Important conventions

- Paths in `input.nmlt` are interpreted relative to the process working directory, normally the repository root.
- The code uses double precision (`kind(1.0d0)`) throughout the shared numerical modules.
- The normalized radial coordinate is internally called `radcor`/`x`; the grid is generated from normalized poloidal flux limits and then mapped to the code’s radial coordinate.
- The magnetic-coordinate arrays use `(poloidal, radial)` indexing in many routines, while perturbation fields are commonly stored as `(toroidal, radial, left/right z-cell boundary)`.
- `mtor` is the toroidal grid count and must be even for the FFT implementation used by the code. `mpol` is the equilibrium poloidal count and is expected to be odd.
- A “tube” is a toroidal/MPI decomposition group. `ntube` is not merely a physics parameter: it participates directly in communicator construction and process-count checks.
