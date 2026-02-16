#include <omp.h>
#include "parallelwork.h"

void work_alloc(parallel_work_t& work, int maxN, int sendrecvnum, int bufmax, int maxsendcells)
{
    #pragma omp parallel default(none) shared(work)
    #pragma omp single
    work.nthreads = omp_get_num_threads();    
    work.atomfx = rarray<double,2>(work.nthreads, maxN);
    work.atomfy = rarray<double,2>(work.nthreads, maxN);
    work.atomfz = rarray<double,2>(work.nthreads, maxN);
    work.send_buffer_atoms = rarray<atom_t,2>(sendrecvnum, bufmax);
    work.recv_buffer_atoms = rarray<atom_t,2>(sendrecvnum, bufmax);
    work.blocklens = rarray<int,2>(sendrecvnum, maxsendcells);
    work.blockinit = rarray<int,2>(sendrecvnum, maxsendcells);
}

void work_free(parallel_work_t& work)
{
    work.atomfx.clear();
    work.atomfy.clear();
    work.atomfz.clear();
    work.send_buffer_atoms.clear();
    work.recv_buffer_atoms.clear();
    work.blocklens.clear();
    work.blockinit.clear();
    work.nthreads = 0;
}
