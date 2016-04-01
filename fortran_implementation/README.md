# QO_fort
### A newer Fortran implementation of the recursive calculation of Quantum Oscillations

This is also written in an object oriented manner, but now in the much faster (compiled)
Fortran language. Usage is virtually identical to the python implementation. (see main_recursive.f90)

### Description copied from python implementation
An object oriented version of the recursive Green's function approach to calculation Quantum Oscillations. Here we include a 6 band Hamiltonian representing the three t2g orbitals ($d_{xz}$, $d_{yx}$ and $d_{xy}$), and including spin orbit coupling.

Details of the recursive algorithm are given here:
http://diku.dk/forskning/Publikationer/tekniske_rapporter/2010/Inversion-of-Block-Tridiagonal-Matrices.pdf

### General format of the code:
Initialize an object which grabs input parameters from the commnad line
Initialize a Hamiltonian object which has appropriate (tight binding) parameters.
Pass this Hamiltonian object to the Solver object, which implements the recursive algorithm.

### Example usage from command line:
./qo length=16394 mu=-0.5 nsteps=1000 spinOrbit=0.01 minInvField=10 maxInvField=1000