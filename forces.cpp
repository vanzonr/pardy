/* force.cpp */
#include "forces.h"
#include "debug.h"
#include <omp.h>
#include <cassert>
#include <cmath>

void interaction_pairs_alloc(interaction_pairs_t& p, long long Npairsmax)
{
    p.npairs = Npairsmax;
    p.pairi  = rvector<int>(Npairsmax);
    p.pairj  = rvector<int>(Npairsmax);
    p.dj     = rvector<packed_int>(Npairsmax);
}

void interaction_pairs_free(interaction_pairs_t& p)
{
    p.pairi.clear();
    p.pairj.clear();
    p.dj.clear();
}

double computeForces(int N, rvector<atom_t>& atoms, interaction_pairs_t& p, double L, bool accum, parallel_work_t& work)
/* accum: true if the forces need to be added too, false if the forces should be replaced 
   returns the potential enenty.*/
{
    #ifndef NDEBUG
    int nth = 0;
    #pragma omp parallel reduction(max:nth)
    nth = omp_get_num_threads();    
    assert_ge(work.atomfx.extent(0), nth);
    assert_ge(work.atomfy.extent(0), nth);
    assert_ge(work.atomfz.extent(0), nth);
    assert_ge(work.atomfx.extent(1), N);
    assert_ge(work.atomfy.extent(1), N);
    assert_ge(work.atomfz.extent(1), N);
    #endif
    /* initialize energy and forces to zero */
    double Usum = 0;
    if (!accum) {
        #pragma omp parallel for default(none) shared(N,atoms)
        for (int i = 0; i < N; ++i) 
            atoms[i].fx = atoms[i].fy = atoms[i].fz = 0.0;
    } 
    /* determine interactions */
    #pragma omp parallel reduction(+:Usum)
    {
        int c = omp_get_thread_num();
        for (long long k = 0; k < N; ++k)
            work.atomfx[c][k] = work.atomfy[c][k] = work.atomfz[c][k] = 0.0;
        #pragma omp for 
        for (long long k = 0; k < p.npairs; ++k) {
            int i = p.pairi[k];
            int j = p.pairj[k];
            int d = p.dj[k];
            /* determine distance in periodic geometry */
            double dx = atoms[i].rx - atoms[j].rx - L*pack2x(d);
            double dy = atoms[i].ry - atoms[j].ry - L*pack2y(d); 
            double dz = atoms[i].rz - atoms[j].rz - L*pack2z(d);
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < rc*rc) {
                double r2i = 1/r2;
                double r6i = r2i*r2i*r2i;
                double fij = 48*r2i*r6i*(r6i-0.5); /* force multiplier */
                double eij = 4*r6i*(r6i-1); /* potential energy between i and j */
                /* within smooth cutoff region? */
                if (r2 > rcp*rcp) {
                    /* in out part of potential, between rcp and rc,      */
                    /* replace phi by alpha * phi for smoother potential  */
                    /* based on a slight rewriting of the alpha factor to */
                    /*  alpha=1/2-1/4x(x^2-3)                             */
                    /* where                                              */
                    /*  x=(2r-rcp-rc)/(rcp-rc)                            */          
                    double r = sqrt(r2);
                    double x = (2*r-rcp-rc)/(rcp-rc);
                    double alpha  = 0.5-0.25*x*(x*x-3);
                    double dalpha = 1.5*(x*x-1)/(r*(rcp-rc));
                    fij = alpha*fij+dalpha*eij;
                    eij = alpha*eij;
                }
                Usum += eij;
                work.atomfx[c][i] += fij*dx;
                work.atomfy[c][i] += fij*dy;
                work.atomfz[c][i] += fij*dz;
                work.atomfx[c][j] -= fij*dx;
                work.atomfy[c][j] -= fij*dy;
                work.atomfz[c][j] -= fij*dz;
            }
        }
        /* must do reduction of c arrays by hand in openmp, at least before version 4.5 (gcc 6.1?) */
        #pragma omp for
        for (int k = 0; k < N; ++k)
            for (int c = 0; c < work.nthreads; ++c) {
                atoms[k].fx += work.atomfx[c][k];
                atoms[k].fy += work.atomfy[c][k];
                atoms[k].fz += work.atomfz[c][k];
            }   
    }/*end omp parallel*/
    return Usum;
}
