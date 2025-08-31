TEK performs particle in cell simulation of low-frequency electromagnetic turbulence in tokamak plasmas using gyrokinetic model for ions, and drift-kinetic for electrons.
Documents about TEK: https://youjunhu.github.io/research_notes/nonlinear_gyrokinetic_equation.pdf

Libraries requried by TEK: Lapack and FFTW

Comile and run TEK:

In linux terminal, compile TEK:

$ make

Run TEK:

mpirun -n 256  ./TEK


Let me know if you have interst in using TEK.

Although TEK is still a one-person code, it has all the capabilities one expect a modern gyrokinetic code to have: kinetic electrons and electromagnetic perturbations, being capable to simulate kinetic ballooning modes and trapped electron modes.