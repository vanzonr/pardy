///
/// @file atom.cpp
///
/// @brief Implements the atom_t type and its MPI type.
///
/// @author Ramses van Zon
/// @date 2026
///
#include <mpi.h>
#include "atom.h"
#include <array>

MPI_Datatype MPI_ATOM = MPI_BYTE;
static bool MPI_ATOM_defined = false;

MPI_Datatype define_MPI_ATOM()
{
    if (! MPI_ATOM_defined) {
        constexpr int num_atom_properties = 14;
        std::array<int,num_atom_properties> array_of_ones {
            1, 1, 1,
            1, 1, 1,
            1, 1, 1,
            1,
            1,1,1,
            1};
        std::array<MPI_Aint,num_atom_properties> array_of_offsets {
            offsetof(atom_t,rx), offsetof(atom_t,ry), offsetof(atom_t,rz),
            offsetof(atom_t,px), offsetof(atom_t,py), offsetof(atom_t,pz),
            offsetof(atom_t,fx), offsetof(atom_t,fy), offsetof(atom_t,fz),
            offsetof(atom_t,index),
            offsetof(atom_t,cx),offsetof(atom_t,cy),offsetof(atom_t,cz),
            offsetof(atom_t,c)};
        std::array<MPI_Datatype,num_atom_properties> array_of_types {
            MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
            MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
            MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
            MPI_LONG_LONG,
            MPI_INT, MPI_INT, MPI_INT,
            MPI_INT};
        MPI_Datatype MPI_ATOM_INNER;
        MPI_Type_create_struct(num_atom_properties,
                               array_of_ones.data(),
                               array_of_offsets.data(),
                               array_of_types.data(),
                               &MPI_ATOM_INNER);
        MPI_Type_create_resized(MPI_ATOM_INNER, 0, sizeof(atom_t), &MPI_ATOM);
        MPI_Type_commit(&MPI_ATOM);
        MPI_ATOM_defined = true;
    }
    return MPI_ATOM;
    
}
