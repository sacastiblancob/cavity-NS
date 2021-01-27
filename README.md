# 2D NAVIER STOKES
> Bidimensional Solver of Navier Stokes equations with Finite Differences - FORTRAN

2D Navier Stokes finite differences solver with CSC storage and solvers for matrices.

With Fortran - GNU-Fortran compiler (and Matlab originals).

Most relevant things in this solver are the Compressed Sparse Column (CSC) tools developed.

![plot](./images/cavity.png)

## Installation

OS X & GNU-Linux:

```sh
git clone https://github.com/sacastiblancob/cavity_NS.git
```

```sh
make clean; make
```

## Usage example

To run the solver you should have a Fotran compiler (if possible use GNU-Fortran!!), change in the Makefile the lines related with the compiler and user specific configuration (the actual version works perfectly with Gfortran).

Once the Makefile is modified (if you need to), you should compile it by typing the make clean and make commands. The .out executable file will be located in the ./bin/ folder.

Modify the configuration in the ./nsconf.nml file, if you want to change input variables, modifiable ones are explained within this file.

Once you have configured your scenario and compile the solver, you can run it by the following command:

```sh
./bin/ns_df.out
```

Wait to reach the Final Time and all the information you need to know will be printed in the terminal, it will show you values such as the tolerance and number of iterations used by the SOR solver in the pressure step, the number of iterations of Conjugate Gradient Solver in diffusion step, Reynolds number, etc.

The results will be storaged in the ./res/ forder, in the .dat files. You can load them into Paraview by using the "TecPlot Reader" option.

Printed Variables: Coordinates, Velocity in X direction (U), Velocity in Y direction (V), Pressure (P) and magnitude of velocity (U_MAG)

## Folder Contents

./bin/ --> Executable file

./liquid/ --> Boundary and initial condition files (Developing stage)

./matlab/ --> Matlab Files (main = Solver_PF.m)

./mod/ --> Fortran modules

./obj/ --> Fortran objects

./res/ --> Solver results(.dat files)

./src/ --> Fotran source codes

./nsconf.nml --> Configuration file (change input variables here)

./Makefile --> Make compiler file

./README.md --> You are standing here

## Meta
Sergio A. Castiblanco-Ballesteros

Bogota - Colombia

Mail 1: sacastiblancob@unal.edu.co

Mail 2: sergio.castiblanco@javeriana.edu.co

> Free Distribution and Open Source (As it should be!!)


