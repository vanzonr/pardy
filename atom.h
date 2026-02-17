/* atom.h - defines the atom_t type and its mpi type */
#ifndef ATOMH
#define ATOMH

#include <mpi.h>

struct atom_t {
    double    rx, ry, rz; /**< position */
    double    px, py, pz; /**< momentum */
    double    fx, fy, fz; /**< force */
    long long index;      /**< mostly for debugging: which particle are you? */
    int       cx, cy, cz; /**< cell indices in each directorion */
    int       c;          /**< super-cell-index (mapping from c{x,y,z} varies)*/
};

/**
 * MPI data type to exchange atoms between processes
 */
extern MPI_Datatype MPI_ATOM;


/**
 * Function to create and commit the MPI data type MPI_ATOM, which
 * must be done once before use. Call this function just after your
 * MPI_Init() call.
 */
MPI_Datatype define_MPI_ATOM();

#endif
