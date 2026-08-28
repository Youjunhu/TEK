TEK is a delta-f gyrokinetic PIC code, which can simulate low-frequency electromagnetic turbulence in tokamak plasmas using the gyrokinetic model for all species (with electrons in the zero Larmor radius limit).

External libraries needed: Lapack and FFTW

Comile TEK:

$ make compile_one_file

Run TEK:

mpirun -n 512  ./TEK

Input to TEK: a fortran namelist `inmput.nmlt`, which specifies various numerical parameters (e.g. grid point number), filenames of equilibrium magnetic configuration (g-eqdsk file), and species radial profiles (e.g. number density and temperature).


TEK has been benchmarked with the GENE code in the DIII-D cyclone base case for the ITG-KBM transition and ITG-TEM transition.
TEK was also benchmarked with the ORB5 code and TRIMEG-GKX code in the ITPA-EP TAE benchmarking.

Documents about TEK: https://arxiv.org/abs/2608.06764

More details can be found in my personal notes: https://youjunhu.github.io/research_notes/nonlinear_gyrokinetic_equation.pdf

