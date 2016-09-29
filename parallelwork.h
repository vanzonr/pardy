#ifndef PARALLELWORKH
#define PARALLELWORKH

#include "atom.h"

typedef struct parallel_work_t {
    double** atomfx;
    double** atomfy;
    double** atomfz;
    int      nthreads;
    atom_t** send_buffer_atoms;
    atom_t** recv_buffer_atoms;
    int**    blocklens;
    int**    blockinit;
} parallel_work_t;

void work_alloc(parallel_work_t* work, int maxN, int sendrecvnum, int bufmax, int maxsendcells);

void work_free(parallel_work_t* work);

#endif
