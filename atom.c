/* atom.c - implements the atom_t type and its mpi type */
#include <stdbool.h>
#include <mpi.h>
#include "atom.h"

MPI_Datatype MPI_ATOM = MPI_BYTE;

static bool MPI_ATOM_defined = false;

MPI_Datatype define_MPI_ATOM()
{
    if (! MPI_ATOM_defined) {
        MPI_Type_struct(16, 
                       (int[16]) {
                           1,
                           1, 1, 1,
                           1, 1, 1,
                           1, 1, 1,
                           1,
                           1,1,1,
                           1,
                           1},
                       (MPI_Aint[16]) {
                           0,
                           offsetof(atom_t,rx), offsetof(atom_t,ry), offsetof(atom_t,rz),
                           offsetof(atom_t,px), offsetof(atom_t,py), offsetof(atom_t,pz),
                           offsetof(atom_t,fx), offsetof(atom_t,fy), offsetof(atom_t,fz),
                           offsetof(atom_t,index),
                           offsetof(atom_t,cx),offsetof(atom_t,cy),offsetof(atom_t,cz),
                           offsetof(atom_t,c),
                           sizeof(atom_t)},
                       (MPI_Datatype[16]) {
                           MPI_LB,
                           MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                           MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                           MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
                           MPI_LONG_LONG,
                           MPI_INT, MPI_INT, MPI_INT,
                           MPI_INT,
                           MPI_UB},
                       &MPI_ATOM);
        //atom_t testatoms[2];
        //MPI_Type_contiguous(&(testatoms[1])-&(testatoms[0]), MPI_BYTE, &MPI_ATOM);
        MPI_Type_commit(&MPI_ATOM);
        MPI_ATOM_defined = true;
    }
    return MPI_ATOM;
    
}
