#include "mpicommunication.h"
#include <cassert>
#include <algorithm>
#include <ndmalloc.h>
#include "cells.h"
#include "debug.h"
#include "global.h"
#include "estimates.h"

static double periodic(double u, double L)
{
    while (u < -0.5*L)
        u += L;
    while (u >= 0.5*L)
        u -= L;
    return u;
} 

#define include13(x)            (((x)>=13)?(x)-1:(x))
#define skip13(x)               (((x)>=13)?(x)+1:(x))
#define SEND_CELL_TO(direction) do{int so=include13(direction);if(!sys.neighborsendblock[direction]){assert_lt(blockcount[so],maxsendcells);work.blockinit[so][blockcount[so]]=sys.start_of_cell[c];work.blocklens[so][blockcount[so]]=sys.n_in_cell[c];blockcount[so]++;}}while(false)

void detect_and_send_ghost_particles(rvector<atom_t>& atoms, size_t bufmax,
                                     atom_t* recv_buffer_atoms[],
                                     MPI_Request send_request[],
                                     MPI_Request recv_request[],
                                     int& nisends,
                                     int& nirecvs,
                                     system_t& sys,
                                     parallel_work_t& work)
{
    /* build index lists */
    #ifndef NDEBUG
    int maxsendcells = std::max({sys.nc[0]*sys.nc[1], sys.nc[0]*sys.nc[2], sys.nc[1]*sys.nc[2]});
    #endif
    int blockcount[SENDNUM] = {0};
    assert(ndsize(work.blocklens,0) >= SENDNUM);
    assert(ndsize(work.blockinit,0) >= SENDNUM);
    assert(ndsize(work.blocklens,1) >= maxsendcells);
    assert(ndsize(work.blockinit,1) >= maxsendcells);
    assert(ndsize(recv_buffer_atoms,0) >= RECVNUM);
    assert(ndsize(recv_buffer_atoms,1) >= bufmax);
    for (int c = 0; c < sys.ncprod; ++c) {
        int cx = super_cell_index_2_cx(c, sys.minc, sys.nc);
        int cy = super_cell_index_2_cy(c, sys.minc, sys.nc);
        int cz = super_cell_index_2_cz(c, sys.minc, sys.nc);
        /* send corner cells: */
        if (cx == sys.minc[0] && cy == sys.minc[1] && cz == sys.minc[2]) 
            SEND_CELL_TO(BACKBOTTOMLEFT);
        if (cx == sys.maxc[0] && cy == sys.minc[1] && cz == sys.minc[2])
            SEND_CELL_TO(FRONTBOTTOMLEFT);
        if (cx == sys.minc[0] && cy == sys.maxc[1] && cz == sys.minc[2])
            SEND_CELL_TO(BACKBOTTOMRIGHT);
        if (cx == sys.maxc[0] && cy == sys.maxc[1] && cz == sys.minc[2])
            SEND_CELL_TO(FRONTBOTTOMRIGHT);
        if (cx == sys.minc[0] && cy == sys.minc[1] && cz == sys.maxc[2]) 
            SEND_CELL_TO(BACKTOPLEFT);
        if (cx == sys.maxc[0] && cy == sys.minc[1] && cz == sys.maxc[2])
            SEND_CELL_TO(FRONTTOPLEFT);
        if (cx == sys.minc[0] && cy == sys.maxc[1] && cz == sys.maxc[2])
            SEND_CELL_TO(BACKTOPRIGHT);
        if (cx == sys.maxc[0] && cy == sys.maxc[1] && cz == sys.maxc[2])
            SEND_CELL_TO(FRONTTOPRIGHT);
        /* send edge cells: */
        if (cy == sys.minc[1] && cz == sys.minc[2])
            SEND_CELL_TO(LEFTBOTTOM);
        if (cy == sys.maxc[1] && cz == sys.minc[2])
            SEND_CELL_TO(RIGHTBOTTOM);
        if (cy == sys.minc[1] && cz == sys.maxc[2])
            SEND_CELL_TO(LEFTTOP);
        if (cy == sys.maxc[1] && cz == sys.maxc[2])
            SEND_CELL_TO(RIGHTTOP);
        if(cx == sys.minc[0] && cz == sys.minc[2])
            SEND_CELL_TO(BACKBOTTOM);
        if (cx == sys.maxc[0] && cz == sys.minc[2])
            SEND_CELL_TO(FRONTBOTTOM);
        if (cx == sys.minc[0] && cz == sys.maxc[2])
            SEND_CELL_TO(BACKTOP);
        if (cx == sys.maxc[0] && cz == sys.maxc[2])
            SEND_CELL_TO(FRONTTOP);
        if (cx == sys.minc[0] && cy == sys.minc[1])
            SEND_CELL_TO(BACKLEFT);
        if (cx == sys.maxc[0] && cy == sys.minc[1])
            SEND_CELL_TO(FRONTLEFT);
        if (cx == sys.minc[0] && cy == sys.maxc[1])
            SEND_CELL_TO(BACKRIGHT);
        if (cx == sys.maxc[0] && cy == sys.maxc[1])
            SEND_CELL_TO(FRONTRIGHT);
        /* send face cells: */
        if (cx == sys.minc[0])
            SEND_CELL_TO(BACK);
        if (cx == sys.maxc[0])
            SEND_CELL_TO(FRONT);
        if (cy == sys.minc[1])
            SEND_CELL_TO(LEFT);
        if (cy == sys.maxc[1])
            SEND_CELL_TO(RIGHT);
        if (cz == sys.minc[2])
            SEND_CELL_TO(BOTTOM);
        if (cz == sys.maxc[2])
            SEND_CELL_TO(TOP);
    }
    /* use these to send using mpi_datatype */
    MPI_Datatype indexedsend[SENDNUM];
    nisends = 0;    
    for (int sendorder = 0; sendorder < SENDNUM; ++sendorder) {
        int target = skip13(sendorder);
        if (sys.neighborrank[target] != sys.rank && ! sys.neighborsendblock[target]) {
            #ifndef NDEBUG
            size_t nsend = 0;
            for (int i = 0; i < blockcount[sendorder]; ++i)
               nsend += work.blocklens[sendorder][i];
            assert_le(nsend, bufmax);
            #endif
            MPI_Type_indexed(blockcount[sendorder], work.blocklens[sendorder], work.blockinit[sendorder], MPI_ATOM, &(indexedsend[nisends]));
            MPI_Type_commit(&(indexedsend[nisends]));
            MPI_Isend(atoms.data(), 1, indexedsend[nisends], sys.neighborrank[target], sendorder, sys.comm, &(send_request[nisends]));
            nisends += 1;
        } else 
            assert(blockcount[sendorder] == 0);
    }
    /* also setup receives */
    nirecvs = 0;
    for (int recvorder = 0; recvorder < RECVNUM; ++recvorder) {
        int source = skip13(RECVNUM-1-recvorder);
        if (sys.neighborrank[source] != sys.rank
            && ! sys.neighborsendblock[26-source] ) { /* receive 'source' is blocked if send '26-source' is blocked */
            MPI_Irecv(recv_buffer_atoms[nirecvs], bufmax, MPI_ATOM, sys.neighborrank[source], MPI_ANY_TAG, sys.comm, &(recv_request[nirecvs]));
            nirecvs += 1;
        }
    }
    for (int i = 0; i < nisends; ++i) 
        MPI_Type_free(&(indexedsend[i])); 
}

int wait_for_particles(int num_send_requests, MPI_Request* send_requests,
                       int num_recv_requests, MPI_Request* recv_requests,
                       atom_t* recv_buffer_atoms[], rvector<atom_t> atoms_new)
{
    MPI_Status  status[num_recv_requests+num_send_requests];
    MPI_Request all_requests[num_send_requests+num_recv_requests];
    for (int i = 0; i < num_recv_requests; ++i)
        all_requests[i] = recv_requests[i];
    for (int i = 0; i < num_send_requests; ++i)
        all_requests[num_recv_requests+i] = send_requests[i];    
    MPI_Waitall(num_recv_requests+num_send_requests, all_requests, status);
    size_t count_new = 0;
    for (int i = 0; i < num_recv_requests; ++i) {
        int count = 0;
        MPI_Get_count(&(status[i]), MPI_ATOM, &count);
        #ifndef NDEBUG
        size_t bufmax = ndsize(recv_buffer_atoms,1);
        assert_le(count, bufmax);
        #endif
        for (int j = 0; j < count; ++j) {
            atoms_new[count_new+j].index = recv_buffer_atoms[i][j].index;
            atoms_new[count_new+j].rx = recv_buffer_atoms[i][j].rx;
            atoms_new[count_new+j].ry = recv_buffer_atoms[i][j].ry;
            atoms_new[count_new+j].rz = recv_buffer_atoms[i][j].rz;
            atoms_new[count_new+j].px = recv_buffer_atoms[i][j].px;
            atoms_new[count_new+j].py = recv_buffer_atoms[i][j].py;
            atoms_new[count_new+j].pz = recv_buffer_atoms[i][j].pz;
            atoms_new[count_new+j].fx = recv_buffer_atoms[i][j].fx;
            atoms_new[count_new+j].fy = recv_buffer_atoms[i][j].fy;
            atoms_new[count_new+j].fz = recv_buffer_atoms[i][j].fz;
        }
        count_new += count;
    }
    return count_new;
}

void detect_and_send_exchange_particles(rvector<atom_t>& atoms, size_t bufmax,
                                        atom_t* send_buffer_atoms[],
                                        atom_t* recv_buffer_atoms[],
                                        MPI_Request send_requests[],
                                        MPI_Request recv_requests[],
                                        int& nisends,
                                        int& nirecvs,
                                        system_t& sys)
{
    int bufcount[SENDNUM] = {0};
    #pragma omp parallel for
    for (int i = 0; i < sys.N; ++i) {
        int ix = 0;
        int iy = 0;
        int iz = 0;
        if (sys.np[0] > 1) {
            if (atoms[i].rx < sys.origin[0])
                ix = -1;
            else if (atoms[i].rx >= sys.origin[0] + sys.localL[0])
                ix = +1;
        }
        if (sys.np[1] > 1) {
            if (atoms[i].ry < sys.origin[1])
                iy = -1;
            else if (atoms[i].ry >= sys.origin[1] + sys.localL[1])
                iy = +1;
        }
        if (sys.np[2] > 1) {
            if (atoms[i].rz < sys.origin[2])
                iz = -1;
            else if (atoms[i].rz >= sys.origin[2] + sys.localL[2])
                iz = +1;
        }
        atoms[i].c = CENTER+ix+3*iy+9*iz; /* temporarily reuse c for new proc where to send */
        /* apply periodic boundary conditions: */
        atoms[i].rx = periodic(atoms[i].rx, sys.L);
        atoms[i].ry = periodic(atoms[i].ry, sys.L);
        atoms[i].rz = periodic(atoms[i].rz, sys.L);
    }
    /* take out atoms */
    for (int i = 0; i < sys.N; ++i) 
        while (atoms[i].c != CENTER && i < sys.N) {
            int sendorder = include13(atoms[i].c);
            send_buffer_atoms[sendorder][bufcount[sendorder]] = atoms[i];
            atoms[i] = atoms[sys.N-1];
            bufcount[sendorder]++;
            sys.N --;
        }
    #ifndef NDEBUG
    for (int c = 0; c < SENDNUM; ++c)  {
        int thisbufcount = bufcount[include13(c)];
        assert_lt(thisbufcount, bufmax);
    }
    #endif
    /* setup sends */
    nisends = 0;
    for (int sendorder = 0; sendorder < SENDNUM; ++sendorder) {
        int target = skip13(sendorder);
        if (sys.neighborrank[target] != sys.rank) {
            MPI_Isend(&(send_buffer_atoms[sendorder][0]), bufcount[sendorder], MPI_ATOM, sys.neighborrank[target], sendorder, sys.comm, &(send_requests[nisends]));
            nisends += 1;
        } else {
            /* omit sends to self */
            assert(bufcount[sendorder] == 0); /* check that the program wasn't going to send anything anyway. */
        }
    }
    /* also setup receives */
    nirecvs = 0;
    for (int recvorder = 0; recvorder < RECVNUM; ++recvorder) {
        int source = skip13(RECVNUM-1-recvorder);
        if (sys.neighborrank[source] != sys.rank) {
            MPI_Irecv(recv_buffer_atoms[nirecvs], bufmax, MPI_ATOM, sys.neighborrank[source], MPI_ANY_TAG, sys.comm, &(recv_requests[nirecvs]));
            nirecvs += 1;
        } else {
            /* omit recv from self */
            send_requests[recvorder] = MPI_REQUEST_NULL;
        }
    }
}

void exchangeParticles(rvector<atom_t>& atoms, system_t& sys, parallel_work_t& work)
{
    MPI_Request  send_requests[SENDNUM];
    MPI_Request  recv_requests[RECVNUM];
    int          nisends = 0;
    int          nirecvs = 0;
    int bufmax = estimateMaxNCrossThroughFace(sys.rho, sys.localL, sys.T0, sys.dt, 0.5*rcp);
    assert(ndsize(work.send_buffer_atoms,0) >= SENDNUM);
    assert(ndsize(work.send_buffer_atoms,1) >= bufmax);
    assert(ndsize(work.recv_buffer_atoms,0) >= RECVNUM);
    assert(ndsize(work.recv_buffer_atoms,1) >= bufmax);
    detect_and_send_exchange_particles(atoms, bufmax,
                                       work.send_buffer_atoms, work.recv_buffer_atoms,
                                       send_requests, recv_requests,
                                       nisends, nirecvs, sys);
    size_t exchange_count = wait_for_particles(nisends, send_requests,
                                               nirecvs, recv_requests,
                                               work.recv_buffer_atoms,
                                               atoms.slice(sys.N, atoms.size()));
    sys.N += exchange_count;
}

/* integration and measurement */
void distribute(system_t& sys)
{
    /* divide processors first */
    int npmax = (int)(sys.L/rc);
    float FnpUse = 0;
    float z = 1.0f/(npmax*npmax+npmax*npmax+npmax*npmax);
    /* brute force way of getting the most uniform 3d distribution for
       processor topology */
    for (int iz = 1; iz <= npmax; ++iz) {       
        for (int iy = 1; iy <= npmax; ++iy) {
            for (int ix = 1; ix <= npmax; ++ix) {
                float tryFnpUse = ix*iy*iz
                    -z*((ix-iy)*(ix-iy)+(iz-iy)*(iz-iy)+(ix-iz)*(ix-iz));
                if (tryFnpUse <= sys.nprocs && tryFnpUse > FnpUse) {
                    FnpUse = tryFnpUse;
                    sys.np[0] = ix;
                    sys.np[1] = iy;
                    sys.np[2] = iz;
                }
            }
        }
    }
    /* sanity check */
    if (sys.np[0]*sys.np[1]*sys.np[2] != sys.nprocs) {
        if (i_am_root) {
            fprintf(stderr, "\nPARDY PARALLEL CONSISTENCY ERROR: Cannot use %d processors or each domain would be smaller than the interaction range (%d processors would work).\n\nPARDY PARALLEL PARAMETER ERROR: Reduce the number of processors!\n\n", sys.nprocs, sys.np[0]*sys.np[1]*sys.np[2]);
        }
        MPI_Abort(sys.comm, 6);
    }
    /* make a new communicator with the right topology */
    int periods[DIM] = {1,1,1};
    MPI_Comm newcomm;
    MPI_Cart_create(sys.comm, DIM, sys.np, periods, true, &newcomm);
    /* reestablish rank */
    sys.comm = newcomm;
    MPI_Comm_rank(sys.comm, &(sys.rank));
    if (sys.rank != global_rank) 
        printf("# Warning: In create_cart, rank was changed from %d to %d\n",
               global_rank, sys.rank);
    if (i_am_root) 
        printf("# (npx npy npz) = (%d %d %d)\n", sys.np[0], sys.np[1], sys.np[2]);
    /* determine neighboring domains */
    int me[DIM];
    MPI_Cart_coords(sys.comm, sys.rank, DIM, me);    
    for (int d = 0; d < DIM; ++d) {
        sys.localL[d]   = sys.L/sys.np[d];
        sys.origin[d]   = -0.5*sys.L + me[d]*sys.localL[d];
        if (sys.usecells)
            sys.nc[d] = (int)(sys.localL[d]/rc);
        else
            sys.nc[d] = 1;
        sys.totnc[d]    = sys.nc[d] * sys.np[d];
        sys.minc[d]     = me[d] * sys.nc[d];
        sys.maxc[d]     = sys.minc[d] + sys.nc[d] - 1;
        sys.cellsize[d] = sys.localL[d]/sys.nc[d];
    }
    for (int ix = -1; ix <= 1; ++ix)
        for (int iy = -1; iy <= 1; ++iy)
            for (int iz = -1; iz <= 1; ++iz) {
                int iplace[DIM] = {me[0]+ix,me[1]+iy,me[2]+iz};
                MPI_Cart_rank(sys.comm,
                              iplace,
                              &(sys.neighborrank[CENTER+ix+3*iy+9*iz]));
                bool block = false;
                if (sys.np[0] == 1 && ix != 0)
                    block = true;
                if (sys.np[1] == 1 && iy != 0)
                    block = true;
                if (sys.np[2] == 1 && iz != 0)
                    block = true;
                if (sys.np[0] == 2 && sys.nc[0] == 1 && ix < 0)
                    block = true;
                if (sys.np[1] == 2 && sys.nc[1] == 1 && iy < 0)
                    block = true;
                if (sys.np[2] == 2 && sys.nc[2] == 1 && iz < 0)
                    block = true;
                sys.neighborsendblock[CENTER+ix+3*iy+9*iz] = block;
            }
    sys.neighborsendblock[CENTER] = true;
    sys.ncprod = sys.nc[0] * sys.nc[1] * sys.nc[2];
    sys.start_of_cell = (int*)ndmalloc(sizeof(sys.start_of_cell[0]), 1, sys.ncprod);
    sys.n_in_cell     = (int*)ndmalloc(sizeof(sys.n_in_cell[0]), 1, sys.ncprod);
    if (i_am_root) 
        printf("# (ncx ncy ncz) = (%d %d %d)\n", sys.totnc[0], sys.totnc[1], sys.totnc[2]);
}
