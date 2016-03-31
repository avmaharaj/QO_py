# QO_py
### A python implementation of the recursive calculation of Quantum Oscillations

An object oriented version of the recursive Green's function approach to calculation Quantum Oscillations. Here we include a 6 band Hamiltonian representing the three t2g orbitals ($d_{xz}$, $d_{yx}$ and $d_{xy}$), and including spin orbit coupling.

Details of the recursive algorithm are given here:
http://diku.dk/forskning/Publikationer/tekniske_rapporter/2010/Inversion-of-Block-Tridiagonal-Matrices.pdf

### General format of the code:
Initialize a Hamiltonian object which has appropriate (tight binding) parameters.
Pass this Hamiltonian object to the Solver object, which implements the recursive algorithm.
