#ifndef MPICOMMUNICATIONH
#define MPICOMMUNICATIONH

#include "atom.h"
#include "system.h"
#include "parallelwork.h"

#define SENDNUM 26
#define RECVNUM 26

void detect_and_send_ghost_particles(atom_t* atoms, size_t bufmax,
                                     atom_t* recv_buffer_atoms[],
                                     MPI_Request* send_request,
                                     MPI_Request* recv_request,
                                     int* nisends,
                                     int* nirecvs,
                                     system_t* sys,
                                     parallel_work_t* work);

int wait_for_particles(int num_send_requests, MPI_Request* send_requests,
                       int num_recv_requests, MPI_Request* recv_requests,
                       atom_t* recv_buffer_atoms[], atom_t atoms_new[]);

void detect_and_send_exchange_particles(atom_t atoms[], size_t bufmax,
                                        atom_t* send_buffer_atoms[],
                                        atom_t* recv_buffer_atoms[],
                                        MPI_Request* send_requests,
                                        MPI_Request* recv_requests,
                                        int* nisends,
                                        int* nirecvs,
                                        system_t* sys);

void exchangeParticles(atom_t atoms[], system_t* sys, parallel_work_t* work);

void distribute(system_t* sys);

#endif
