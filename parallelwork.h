///
/// @file parallelwork.h
///
/// @brief TODO.
///
/// @author Ramses van Zon
/// @date 2026
///
#ifndef PARALLELWORKH
#define PARALLELWORKH

#include "atom.h"
#include <rarray>

struct parallel_work_t {

    rarray<double,2> atomfx;
    rarray<double,2> atomfy;
    rarray<double,2> atomfz;
    int              nthreads;
    rarray<atom_t,2> send_buffer_atoms;
    rarray<atom_t,2> recv_buffer_atoms;
    rarray<int,2>    blocklens;
    rarray<int,2>    blockinit;
};

void work_alloc(parallel_work_t& work, int maxN, int sendrecvnum, int bufmax, int maxsendcells);

void work_free(parallel_work_t& work);

#endif
