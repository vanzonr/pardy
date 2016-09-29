/* 
 *  pardy.c
 *
 *  Example of a PARallel molecular DYnamics simulation in C, using
 *  mpi and openmp.
 *
 *  Input (from standard input or filename given on command line): 
 *   N = <number of particles>
 *   rho = <density of the system>
 *   T = <initial temperature (standard deviation of the velocities)>
 *   runtime = <runtime>
 *   dt = <time step>
 *   seed = <random number generator seed>
 *   equil = <equilibration time>
 *   usecells
 *   (when the latter flag is not present, cells will not be used
 *   within processes)
 *
 *  Output (from standard output): lines containing:
 *   time, energy E, pot. en. U, kin. en. K, temperature T, fluctuations, walltime-per-step, walltime-overall
 *
 *  Notes:
 *  - To compile: make
 *  - To run: ./pardy input.ini
 *    where input parameters are listed in the file "input.ini", and
 *    output is sent to stdout.
 *  - Appropriate values for the equilibrium time are best found by
 *    doing a short run and seeing when the potential energy has reach
 *    a stationary value.
 *  - All reported energies values are divided by the number of particles N.
 *  - Fluctuations are the root mean square of E-<E> with <E> the mean energy E.
 *
 *  Ramses van Zon, 13 November 2008
 *  - In Nov 25 2009:
 *    + condensed code
 *    + OpenMP version
 *  - In Aug 2016:
 *    + Renamed to ljhpc, included mpi header (not used yet)
 *  - In Sep 2016:
 *    + Added cells to determine interacting pairs.
 *    + Added MPI
 *    + Renamed to pardy
 *    + Change input format and modularized
 */

#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>
#include <ndmalloc.h>
#include "lcg.h"
#include "lattice.h"
#include "debug.h"
#include "system.h"
#include "atom.h"
#include "inifile.h"
#include "global.h"
#include "estimates.h"
#include "parallelwork.h"
#include "forces.h"
#include "cells.h"
#include "mpicommunication.h"

void sortParticlesAndComputeForces(int *N, atom_t atoms[], interaction_pairs_t* p, system_t* sys, parallel_work_t* work)
{
    /* sort particles into cells */
    cellDivide(atoms, sys);
    /* send ghost particles */
    int         bufmax = estimateMaxNGhostPerFaceNeighbor(sys->rho, sys->localL, sys->cellsize);
    MPI_Request send_requests[SENDNUM];
    MPI_Request recv_requests[RECVNUM];
    int         sendnum = 0;
    int         recvnum = 0;
    assert(ndsize(work->recv_buffer_atoms,0) >= RECVNUM);
    assert(ndsize(work->recv_buffer_atoms,1) >= bufmax);
    detect_and_send_ghost_particles(atoms, bufmax,
                                    work->recv_buffer_atoms,
                                    send_requests,
                                    recv_requests,
                                    &sendnum, &recvnum,
                                    sys, work);
    /* Part of the force calculation involving local particles
       (computation overlaps with communication) */
    findPairsFromCells(p, sys);
    sys->U = computeForces(*N, atoms, p, sys->L, false, work);
    /* Collect neighbors */
    atom_t* ghost_atoms = atoms + *N;
    size_t  ghost_count = wait_for_particles(sendnum, send_requests,
                                             recvnum, recv_requests,
                                             work->recv_buffer_atoms, ghost_atoms);    
    size_t NselfPlusGhosts = *N+ghost_count;
    #ifndef NDEBUG
    size_t atomsCapacity = ndsize(atoms,0);
    assert_le(NselfPlusGhosts, atomsCapacity);
    #endif
    /* Neighbor part of force calculation */
    findGhostPairsFromCells(p, ghost_count, ghost_atoms, *N, sys);
    sys->U += 0.5*computeForces(NselfPlusGhosts, atoms, p, sys->L, true, work);
}

void initialize(atom_t atoms[], interaction_pairs_t* p, system_t* sys, parallel_work_t* work)
{
    double scale;
    int i;
    /* generate positions */
    /* initialize lattice position generator, then pick at most Ntot (will be less if parallel) */
    lattice_t lat = init_partial_lattice(sys->L, sys->Ntot, sys->origin, sys->localL);
    for (i = 0; i < sys->Ntot; ++i) {
        long long index = make_lattice_position(&lat, &(atoms[i].rx), &(atoms[i].ry), &(atoms[i].rz));
        if (index != NO_LATTICE_POSITION) {
            atoms[i].index = index;
        } else {
            break; /* no more points */
        }
    }
    sys->N = i;
    /* report the number of particles held by each process */
    int Nall[global_size];
    MPI_Allgather(&(sys->N), 1, MPI_INT, Nall, 1, MPI_INT, MPI_COMM_WORLD);
    long long checkNtot = 0;
    for (int i = 0; i < global_size; ++i) {
        if (i_am_root)
            printf("# N[%d]=%d\n", i, Nall[i]);
        checkNtot += Nall[i];
    }
    /* sanity check */
    if (checkNtot != sys->Ntot) {
        if (i_am_root)
            fprintf(stderr, "\nPARDY PARALLEL CONSISTENCY ERROR: Combined number of particles in all processes (%lld) does not match parameter N (%lld)\n", checkNtot, sys->Ntot);
        MPI_Abort(sys->comm, 7);
    }
    /* generate momenta */
    lcg_t rng = lcg_init_rand48(sys->seed);  /* initialize random number
                                                generator used in gaussian
                                                instead of "srand48(seed)" */
    scale = sqrt(sys->T0);
    sys->K = 0;
    /* must do some skipping (4 random numbers per particle)*/
    lcg_skip(&rng, 4*atoms[0].index);
    for (int i = 0; i < sys->N; ++i){
        atoms[i].px = scale*lcg_normal(&rng);
        atoms[i].py = scale*lcg_normal(&rng);
        atoms[i].pz = scale*lcg_normal(&rng);
        lcg_normal(&rng); /* must have even number of normal random numbers */
        sys->K += atoms[i].px*atoms[i].px + atoms[i].py*atoms[i].py + atoms[i].pz*atoms[i].pz;
        if (i<(sys->N)-1) {
            int index_diff = atoms[i+1].index - atoms[i].index;
            if (index_diff > 1)
                lcg_skip(&rng, 4*(index_diff-1));
        }
    }
    /* compute instantaneous kinetic temperature and energy */
    sys->K *= 0.5;
    sortParticlesAndComputeForces(&(sys->N), atoms, p, sys, work);
    double Utot, Ktot;
    MPI_Reduce(&(sys->U), &Utot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(sys->K), &Ktot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    /* report results */
    if (i_am_root) {
        printf("#  step    time    E       U       K       T    <[E-<E>]^2>  walltime(s) cum.walltime(s)\n");
        printf("%7d %7.3f %7.3f %7.3f %7.3f %7.3f   .....         0.000      0.000\n",
               0, 0., (Ktot+Utot) / (sys->Ntot), Utot / (sys->Ntot), Ktot / (sys->Ntot), (2*Ktot)/(DIM*(sys->Ntot)));
    }
}

/* Verlet integration step */
void integrateStep(system_t* sys, atom_t atoms[], interaction_pairs_t* p, parallel_work_t* work)
{
    /* first momentum half step (forces assume to be computed in
       initialization or in the last step) */
    #pragma omp parallel for
    for (int i = 0; i < sys->N; ++i) {
        atoms[i].px += 0.5*sys->dt*atoms[i].fx;
        atoms[i].py += 0.5*sys->dt*atoms[i].fy;
        atoms[i].pz += 0.5*sys->dt*atoms[i].fz;
    }
    /* full free motion step */
    #pragma omp parallel for
    for (int i = 0; i<sys->N; ++i) {
        atoms[i].rx = atoms[i].rx + sys->dt*atoms[i].px;
        atoms[i].ry = atoms[i].ry + sys->dt*atoms[i].py;
        atoms[i].rz = atoms[i].rz + sys->dt*atoms[i].pz;
    }
    /* send atoms that have digressed too far to other processes (also enforces pbc)*/
    exchangeParticles(atoms, sys, work);
    /* positions were changed, so recompute the forces */
    sortParticlesAndComputeForces(&(sys->N), atoms, p, sys, work);
    /* final momentum half-step */
    double K2 = 0;
    #pragma omp parallel for reduction(+:K2)
    for (int i = 0; i < sys->N; ++i) {
        atoms[i].px+=0.5*sys->dt*atoms[i].fx;
        atoms[i].py+=0.5*sys->dt*atoms[i].fy;
        atoms[i].pz+=0.5*sys->dt*atoms[i].fz;
        K2 += atoms[i].px*atoms[i].px + atoms[i].py*atoms[i].py + atoms[i].pz*atoms[i].pz;
    }
    /* finish computing K */
    sys->K = 0.5*K2;
}
    
/* integration and measurement */
void run(system_t* sys)
{
    const int Nmax = estimateNmax(sys->N, sys->rho, sys->localL, sys->cellsize);
    const int bufmax_exchange = estimateMaxNCrossThroughFace(sys->rho, sys->localL, sys->T0, sys->dt, 0.5*rcp);
    const int bufmax_ghost = estimateMaxNGhostPerFaceNeighbor(sys->rho, sys->localL, sys->cellsize);
    const int bufmax = bufmax_exchange>bufmax_ghost?bufmax_exchange:bufmax_ghost;
    const int maxsendcells = max3(sys->nc[0]*sys->nc[1], sys->nc[0]*sys->nc[2], sys->nc[1]*sys->nc[2]);
    const long long Npairsmax = estimateNpairsmax(sys->N, sys->rho, sys->localL, sys->cellsize);
    interaction_pairs_t p;
    parallel_work_t work;
    atom_t*   atoms  = ndmalloc(sizeof(atom_t), 1, Nmax);
    int       numSteps = (int)(sys->runtime/sys->dt+0.5);
    int       equilSteps = (int)(sys->equil/sys->dt+0.5);   
    int       count;                 /* counts time steps */
    int       numPoints = 0;         /* counts measurements */
    double    sumE = 0;              /* total energy accumulated over steps */
    double    sumE2 = 0;             /* total energy squared accumulated */
    double    avgE, avgE2, fluctE;   /* average energy, its square, and its fluctuations */
    interaction_pairs_alloc(&p, Npairsmax);
    work_alloc(&work, Nmax, SENDNUM, bufmax, maxsendcells);
    /* draw initial conditions */
    initialize(atoms, &p, sys, &work);
    double  walltimestart = MPI_Wtime();
    double  walltimelast = walltimestart;
    /* perform equilibration steps */
    for (count = 0; count < equilSteps; ++count) {
        integrateStep(sys, atoms, &p, &work);
        double Utot, Ktot;
        MPI_Reduce(&(sys->U), &Utot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&(sys->K), &Ktot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (i_am_root) {
            double t = MPI_Wtime();
            printf("%7d %7.3f %7.3f %7.3f %7.3f %7.3f   .....    %10.3lf %10.3lf\n",
                   count+1, (count+1)*sys->dt, (Ktot+Utot)/sys->Ntot, Utot/sys->Ntot, Ktot/sys->Ntot, (2*Ktot)/(DIM*sys->Ntot), t-walltimelast, t-walltimestart );
            walltimelast = t;
        }
    }    
    /* perform rest of time steps */
    for (count = equilSteps; count<numSteps; ++count) {
        integrateStep(sys, atoms, &p, &work);
        double Utot, Ktot;
        MPI_Reduce(&(sys->U), &Utot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&(sys->K), &Ktot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (i_am_root) {
            double Etot = Ktot + Utot;
            sumE  += Etot;
            sumE2 += Etot * Etot;
            ++numPoints;
            avgE   = sumE/numPoints;
            avgE2  = sumE2/numPoints;
            if (avgE2 > avgE*avgE) {
                fluctE = sqrt(avgE2-avgE*avgE);
            } else {  /* correct for round-off */
                fluctE = 0.0;
            }
            double t = MPI_Wtime();
            printf("%7d %7.3f %7.3f %7.3f %7.3f %7.3f %10.2e %10.3lf %10.3lf\n",
                   count+1, (count+1)*sys->dt, Etot/sys->Ntot, Utot/sys->Ntot, Ktot/sys->Ntot, (2*Ktot)/(DIM*sys->Ntot), fluctE/sys->Ntot, t-walltimelast, t-walltimestart);
            walltimelast = t;
        }
    }
    interaction_pairs_free(&p);
    work_free(&work);
    ndfree(atoms);
    ndfree(sys->n_in_cell);
    ndfree(sys->start_of_cell);
}

/* program start */
int main(int argc, char* argv[])
{
    /* start mpi with thread support */
    int required_support_for_threads = MPI_THREAD_FUNNELED, provided_support_for_threads;
    MPI_Init_thread(&argc, &argv, required_support_for_threads, &provided_support_for_threads);
    assert_le(provided_support_for_threads, required_support_for_threads);
    /* who am i and how many of us are there? */
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    /* new mpi data types */
    define_MPI_ATOM();
    define_MPI_PARAMETERS();
    /* debug? */
    ENTERDEBUGGER;
    /* report parallel setup */
    if (i_am_root) {
        printf("#");
        if (global_size != 1)
            printf("Number of processes = %d, ", global_size);
        /* get number of threads available:  */
        #pragma omp parallel default(none)
        #pragma omp single
        printf("Number of threads = %d\n", omp_get_num_threads());
    }
    /**/
    system_t sys;
    if (i_am_root) {
        FILE* file = (argc>1)?fopen(argv[1],"r"):stdin;
        if (file == NULL) 
            MPI_Abort(MPI_COMM_WORLD, 1);
        key_value_table_t* ini = calloc(1,sizeof(*ini));
        kvt_read(ini, file);
        sys.Ntot     = kvt_lookup_long_long(ini, "N");
        sys.rho      = kvt_lookup_double(ini, "rho");
        sys.T0       = kvt_lookup_double(ini, "T");
        sys.runtime  = kvt_lookup_double(ini, "runtime");
        sys.dt       = kvt_lookup_double(ini, "dt");
        sys.seed     = kvt_lookup_long(ini, "seed");
        sys.equil    = kvt_lookup_double(ini, "equil");
        sys.usecells = kvt_lookup_entry(ini, "usecells") != NULL;
        sys.L        = cbrt(sys.Ntot/sys.rho);
        if (argc > 1 && file != NULL)
            fclose(file);
        free(ini);
        printf("#Ntot=%lld rho=%f L=%f T=%f runtime=%f dt=%f seed=%ld equil=%f usecells=%d\n",
               sys.Ntot,sys.rho,sys.L,sys.T0,sys.runtime,sys.dt,sys.seed,sys.equil,sys.usecells);
    }
    MPI_Bcast(&sys, 1, MPI_PARAMETERS, 0, MPI_COMM_WORLD);
    sys.comm = MPI_COMM_WORLD;
    sys.nprocs = global_size;
    sys.rank = global_rank;
    if (sanity_check(&sys, rc, i_am_root)) {
        distribute(&sys);
        run(&sys);
        MPI_Finalize();
    } else {        
        MPI_Finalize();
        return 1;
    }
    return 0;
}

