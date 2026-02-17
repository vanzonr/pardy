[//]: # \mainpage

PARDY is an example of a PARallel molecular DYnamics simulation in C++, 
using MPI and OpenMP.

# Compiling

To compile, make sure an MPI C++ compiler is available (mpicxx) and type
```
make
```

To generate documentation, type `make doc`.

# Input

Input contains the simulation parameters and are read from standard
input or a filename given on command line:

```
 N = <number of particles>
 rho = <density of the system>
 T = <initial temperature (standard deviation of the velocities)>
 runtime = <runtime>
 dt = <time step>
 seed = <random number generator seed>
 equil = <equilibration time>
 usecells
```
(when the latter flag is not present, cells will not be used
 within processes)

# Output

The current versions of PARDY creates only minimal output. Each time step, it writes a line to console containing:

```
 time, energy E, pot. en. U, kin. en. K, temperature T, fluctuations, walltime-per-step, walltime-overall
```

# Running

To run PARDY, simply type:
```
 mpirun pardy input.ini
```
where input parameters are listed in the file "input.ini", and output is sent to stdout.
    
  - Add "-np N" to the mpirun command to set the number of MPI processes.
  
  - Set OMP_NUM_THREADS to change the number of threads per process.

# Notes

  - Appropriate values for the equilibrium time are best found by
    doing a short run and seeing when the potential energy has reach
    a stationary value.
    
  - All reported energies values are divided by the number of particles N.
  
  - Fluctuations are the root mean square of E-<E> with <E> the mean energy E.

  - This is the C++ version of an older C version. Some parts of the
    code are not quite C++ yet.
