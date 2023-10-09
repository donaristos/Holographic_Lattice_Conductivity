# Holographic Lattice Conductivity

Here you will find the code I wrote to generate the data presented in ["The thermoelectric properties of inhomogeneous holographic lattices"](https://arxiv.org/abs/1409.6875) by A. Donos and J. Gauntlett. The code is in C++ and it required the libraries:
1) [Eigen](https://eigen.tuxfamily.org/index.php?title%253DMain_Page)

   Mostly used to construct and manipulate the Hessian operator in Newton's method implementation

2) [Suitesparse](https://people.engr.tamu.edu/davis/suitesparse.html)

   Used to find the inverse of the Hessian in Newton's method. Especially powerfull in finite differences. 

4) [Mpreal](https://github.com/advanpix/mpreal)

   Multiple precision arithetic library used to perform computations at arbitrary numerical precision.

3) [Boost](https://www.boost.org/)

   Uses the wrapper of the Boost libray to manipulate float128 numbers.