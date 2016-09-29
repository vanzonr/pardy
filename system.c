/* system.h - defines the system_t type and its mpi type */
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include "mpi.h"
#include "system.h"

bool sanity_check(system_t* sys, double rc, bool i_am_root)
{
    /* make sure the system definition makes sense physically and
       computationally. */
    bool sane = true;
    if (sys->Ntot <= 0) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nLJHPC PARAMETER ERROR: Incorrect number of particles (N=%d). "
                    "N must be positive\n",
                    sys->N);
        sane = false;
    }
    if (sys->rho <= 0.0) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nLJHPC PARAMETER ERROR: Incorrect density (rho=%f). "
                    "rho must be positive\n",
                    sys->rho);
        sane = false;
    }
    if (sys->L < 2*rc) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nLJHPC PARAMETER ERROR: Incorrect linear system size (L=%f). "
                    "L must be at least two times the interaction range (rc=%f)\n",
                    sys->L, rc);
        sane = false;
    }
    if (sys->dt <= 0.0) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nLJHPC PARAMETER ERROR: Incorrect time step (dt=%f). "
                    "dt must be positive\n",
                    sys->dt);
        sane = false;
    }
    if (sys->runtime <= 0.0) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nLJHPC PARAMETER ERROR: Incorrect runtime (runtime=%f). "
                    "runtime must be positive\n",
                    sys->runtime);
        sane = false;
    }
    if (sys->seed == 0) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nLJHPC PARAMETER ERROR: Incorrect random number seed (seed=%ld). "
                    "seed must not be zero\n",
                    sys->seed);
        sane = false;
    }
    if  (sys->equil > sys->runtime) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nLJHPC PARAMETER ERROR: Incorrect equilibration time (equil=%f). "
                    "equil must be less then the runtime (%f)\n",
                    sys->equil, sys->runtime);
        sane = false;
    }
    return sane;
}

MPI_Datatype MPI_PARAMETERS = MPI_BYTE;

static bool MPI_PARAMETERS_defined = false;

MPI_Datatype define_MPI_PARAMETERS()
{
    if (! MPI_PARAMETERS_defined) {
        MPI_Type_struct(10,
                        (int[11]) {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
                        (MPI_Aint[11]) {
                            0,
                            offsetof(system_t,Ntot),
                            offsetof(system_t,rho),
                            offsetof(system_t,T0),
                            offsetof(system_t,runtime),
                            offsetof(system_t,dt), 
                            offsetof(system_t,seed),
                            offsetof(system_t,equil),
                            offsetof(system_t,usecells),
                            offsetof(system_t,L),
                            sizeof(system_t)},
                        (MPI_Datatype[11]) {
                            MPI_LB,
                            MPI_LONG_LONG,
                            MPI_DOUBLE,
                            MPI_DOUBLE,
                            MPI_DOUBLE,
                            MPI_DOUBLE,
                            MPI_LONG,
                            MPI_DOUBLE,
                            MPI_C_BOOL,
                            MPI_DOUBLE,
                            MPI_UB},
                       &MPI_PARAMETERS);
        MPI_Type_commit(&MPI_PARAMETERS);
        MPI_PARAMETERS_defined = true;
    }
    return MPI_PARAMETERS;
}

