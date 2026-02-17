/* parallelwork.h */
#ifndef PARALLELWORKH
#define PARALLELWORKH

#include "atom.h"
#include <rarray>

struct parallel_work_t {

    /* double** */ rarray<double,2> atomfx;
    /* double** */ rarray<double,2> atomfy;
    /* double** */ rarray<double,2> atomfz;
    int      nthreads;
    /* atom_t** */ rarray<atom_t,2> send_buffer_atoms;
    /* atom_t** */ rarray<atom_t,2> recv_buffer_atoms;
    /* int**    */ rarray<int,2>    blocklens;
    /* int**    */ rarray<int,2>    blockinit;
};

void work_alloc(parallel_work_t& work, int maxN, int sendrecvnum, int bufmax, int maxsendcells);

void work_free(parallel_work_t& work);

#endif
