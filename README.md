TEK performs particle in cell simulation of low-frequency electromagnetic turbulence in tokamak plasmas using gyrokinetic model for ions, and drift-kinetic for electrons.
Documents about TEK: https://youjunhu.github.io/research_notes/nonlinear_gyrokinetic_equation.pdf

Libraries requried by TEK: Lapack and FFTW

Comile and run TEK:

In linux terminal, compile TEK:

$ make

Run TEK:

mpirun -n 256  ./TEK


TEK has all the capabilities one expects a modern gyrokinetic code to have: kinetic electrons and electromagnetic perturbations, being capable of simulating kinetic ballooning modes and trapped electron modes.

Let me know if you have interest in using TEK.