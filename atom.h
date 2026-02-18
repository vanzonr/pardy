///
/// @file atom.h
///
/// @brief Defines the atom_t type and its MPI type.
///
/// @author Ramses van Zon
/// @date 2026
///
#ifndef ATOMH
#define ATOMH

#include <mpi.h>

struct atom_t {
    double    rx, ry, rz; ///< position
    double    px, py, pz; ///< momentum
    double    fx, fy, fz; ///< force
    long long index;      ///< mostly for debugging: which particle are you?
    int       cx, cy, cz; ///< cell indices in each directorion
    int       c;          ///< super-cell-index (mapping from c{x,y,z} varies)
};

/// @brief MPI data type to exchange atoms between processes
extern MPI_Datatype MPI_ATOM;

/// @brief Create and commit the MPI datatype MPI_ATOM, which must be
///  done once before use. Call this function just after your
///  MPI_Init() call.  Further calls will reuse the existing MPI_ATOM.
///
/// @returns the MPI datatype for atom_t
MPI_Datatype define_MPI_ATOM();

#endif
