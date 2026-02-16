#include <omp.h>
#include <ndmalloc.h>
#include "parallelwork.h"

void work_alloc(parallel_work_t& work, int maxN, int sendrecvnum, int bufmax, int maxsendcells)
{
    #pragma omp parallel default(none) shared(work)
    #pragma omp single
    work.nthreads = omp_get_num_threads();
    work.atomfx = (double**)ndmalloc(sizeof(double), 2, work.nthreads, maxN);
    work.atomfy = (double**)ndmalloc(sizeof(double), 2, work.nthreads, maxN);
    work.atomfz = (double**)ndmalloc(sizeof(double), 2, work.nthreads, maxN);
    work.send_buffer_atoms = (atom_t**)ndmalloc(sizeof(atom_t), 2, sendrecvnum, bufmax);
    work.recv_buffer_atoms = (atom_t**)ndmalloc(sizeof(atom_t), 2, sendrecvnum, bufmax);
    work.blocklens = (int**)ndmalloc(sizeof(int), 2, sendrecvnum, maxsendcells);
    work.blockinit = (int**)ndmalloc(sizeof(int), 2, sendrecvnum, maxsendcells);
}

void work_free(parallel_work_t& work)
{
    ndfree(work.atomfx);
    ndfree(work.atomfy);
    ndfree(work.atomfz);
    ndfree(work.send_buffer_atoms);
    ndfree(work.recv_buffer_atoms);
    ndfree(work.blocklens);
    ndfree(work.blockinit);
    work.atomfx = NULL;
    work.atomfy = NULL;
    work.atomfz = NULL;
    work.nthreads = 0;
    work.send_buffer_atoms = NULL;
    work.recv_buffer_atoms = NULL;
    work.blocklens = NULL;
    work.blockinit = NULL;
}
