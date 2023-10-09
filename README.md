# Holographic Lattice Conductivity

Here you will find the code I wrote to generate the data presented in ["The thermoelectric properties of inhomogeneous holographic lattices"](https://arxiv.org/abs/1409.6875) by A. Donos and J. Gauntlett. The code is in C++ and it required the libraries:
1) [Eigen](https://eigen.tuxfamily.org/index.php?title%253DMain_Page)

   Mostly used to construct and manipulate the Hessian operator in Newton's method implementation. Also used its linear solvers for the flavour of the code that uses high precision numerics.

2) [Intel MKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html)

   Used the functionality of PARDISO to find the inverse of the Hessian in Newton's method when using 64-bit floats. Especially powerfull when the Hessian is sparse (e.g. finite differences).

3) [Mpreal](https://github.com/advanpix/mpreal)

   Multiple precision arithetic library used to perform computations at arbitrary numerical precision.

4) [Boost](https://www.boost.org/)

   Used for higher precision numerics.

To compile the project, I have included a makefile to be used with GNU C++. You might have to modify that depending on where your libraries are located and the distribution of C++ you are using.

The folder "Mathematica" contains the Mathematica notebooks I used to write the equations of motion of Einstein-Maxwell-Dilaton including the DeTurck trick term in Einstein's equations. The subfolder "Background" contains the notebooks concerning the background geometries. The subfolder "Conductivity" contains the notebooks for the equations governing the computation of the optical conductivity.