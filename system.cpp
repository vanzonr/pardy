///
/// @file system.cpp 
///
/// @brief Implementation of system_t functions.
///
/// @author Ramses van Zon
/// @date 2026
///
#include <cstddef>
#include <cstdio>
#include <mpi.h>
#include "system.h"

bool sanity_check(system_t& sys, double rc, bool i_am_root)
{
    // make sure the system definition makes sense physically and computationally.
    bool sane = true;
    if (sys.Ntot <= 0) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nPARDY PARAMETER ERROR: Incorrect number of particles (N=%d). "
                    "N must be positive\n",
                    sys.N);
        sane = false;
    }
    if (sys.rho <= 0.0) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nPARDY PARAMETER ERROR: Incorrect density (rho=%f). "
                    "rho must be positive\n",
                    sys.rho);
        sane = false;
    }
    if (sys.L < 2*rc) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nPARDY PARAMETER ERROR: Incorrect linear system size (L=%f). "
                    "L must be at least two times the interaction range (rc=%f)\n",
                    sys.L, rc);
        sane = false;
    }
    if (sys.dt <= 0.0) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nPARDY PARAMETER ERROR: Incorrect time step (dt=%f). "
                    "dt must be positive\n",
                    sys.dt);
        sane = false;
    }
    if (sys.runtime <= 0.0) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nPARDY PARAMETER ERROR: Incorrect runtime (runtime=%f). "
                    "runtime must be positive\n",
                    sys.runtime);
        sane = false;
    }
    if (sys.seed == 0) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nPARDY PARAMETER ERROR: Incorrect random number seed (seed=%ld). "
                    "seed must not be zero\n",
                    sys.seed);
        sane = false;
    }
    if  (sys.equil > sys.runtime) {
        if (i_am_root) 
            fprintf(stderr,
                    "\nPARDY PARAMETER ERROR: Incorrect equilibration time (equil=%f). "
                    "equil must be less then the runtime (%f)\n",
                    sys.equil, sys.runtime);
        sane = false;
    }
    return sane;
}

MPI_Datatype MPI_PARAMETERS = MPI_BYTE;
static bool MPI_PARAMETERS_defined = false;

MPI_Datatype define_MPI_PARAMETERS()
{
    if (! MPI_PARAMETERS_defined) {
        constexpr int num_parameters = 9;
        std::array<int,num_parameters> array_of_ones {1, 1, 1, 1, 1, 1, 1, 1, 1};
        std::array<MPI_Aint,num_parameters> array_of_offsets {
           offsetof(system_t, Ntot),
           offsetof(system_t, rho),
           offsetof(system_t, T0),
           offsetof(system_t, runtime),
           offsetof(system_t, dt), 
           offsetof(system_t, seed),
           offsetof(system_t, equil),
           offsetof(system_t, usecells),
           offsetof(system_t, L)};
        std::array<MPI_Datatype,num_parameters> array_of_types {
            MPI_LONG_LONG,
            MPI_DOUBLE,
            MPI_DOUBLE,
            MPI_DOUBLE,
            MPI_DOUBLE,
            MPI_LONG,
            MPI_DOUBLE,
            MPI_C_BOOL,
            MPI_DOUBLE};
        MPI_Datatype MPI_PARAMETERS_INNER;
        MPI_Type_create_struct(num_parameters,
                               array_of_ones.data(),
                               array_of_offsets.data(),
                               array_of_types.data(),
                               &MPI_PARAMETERS_INNER);
        MPI_Type_create_resized(MPI_PARAMETERS_INNER, 0, sizeof(system_t), &MPI_PARAMETERS);
        MPI_Type_commit(&MPI_PARAMETERS);
        MPI_PARAMETERS_defined = true;
    }
    return MPI_PARAMETERS;
}
