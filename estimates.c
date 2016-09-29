#include "global.h"
#include "estimates.h"
#include <math.h>
#include "global.h"

int estimateMaxNSelf(int Nmin, double rho, double localL[DIM])
{
 /* Function to estimate how many particle this region will hold at
    most. Based on the volume of the region and its density. Must be
    larger than Nmin.  The formula uses a base of 13 particles (for no
    reason), then allows for 2% systematic variations in the density
    and for 2*sqrt(N) statistical fluctuations in the number of
    particles. */
    double nonfluctuatingN = 13.0 + 1.02*localL[0]*localL[1]*localL[2]*rho;
    return (int)(nonfluctuatingN + 2.0*sqrt(nonfluctuatingN));
}

int estimateMaxNGhostPerFaceNeighbor(double rho, double localL[DIM], double cellsize[DIM])
{
 /* Function to estimate how many particle this region will get from
    an edge neighbor as ghost particles. Based on the maximal surface
    of an edge of the region, the size of a cell, and the average
    system density.  The formula uses a base of 6 particles (for no
    reason), then allows for 55% systematic variations in the density
    and for sqrt(N) statistical fluctuations in the number of
    particles. */
    double nonfluctuatingN = 6.0 + 1.55*max3(localL[0]*localL[1]*cellsize[2],
                                             localL[0]*localL[2]*cellsize[1],
                                             localL[1]*localL[2]*cellsize[0])*rho;
    return (int)(nonfluctuatingN + 2.0*sqrt(nonfluctuatingN));
}

int estimateMaxNCrossThroughFace(double rho, double localL[DIM], double T, double dt, double a)
{
 /* Function to estimate how many particle this region will get from
    an edge neighbor as ghost particles. Based on the temperature, the
    time step, and the range of the potential, as well as the maximal
    surface of an edge of the region.  There is almost no basis to
    this formulat. It uses a base of 1200 particles, scales with
    temperature and time step as a square root, allows with sqrt

    dependency of the result.  a*a is used as an estimate of the collisional cross section. 
    This is inherently the shakiest estimate of the bunch, as initial
    conditions can lead to a large number of exchanged particles
    initially.
*/
    double nonfluctuatingN = 1200 + 1200*sqrt(1+T*dt)/(a*a); 
    return (int)(nonfluctuatingN + 8.0*sqrt(nonfluctuatingN)+5.0*sqrt(rho*localL[0]*localL[1]*localL[2]));
}

int estimateMaxNGhostPerEdgeNeighbor(double rho, double localL[DIM], double cellsize[DIM])
{
 /* Function to estimate how many particle this region will get from a
    neighbor as ghost particles. Based on the maximal length of an
    edge of the region, the size of a cell, and the average system
    density.  The formula uses a base of 6 particles (for no reason),
    then allows for 55% systematic variations in the density and for
    sqrt(N) statistical fluctuations in the number of particles. */
    double nonfluctuatingN = 6.0 + 1.55*max3(
            localL[0]*cellsize[1]*cellsize[2],
            localL[1]*cellsize[0]*cellsize[2],
            localL[2]*cellsize[1]*cellsize[0])*rho;
    return (int)(nonfluctuatingN + 2.0*sqrt(nonfluctuatingN));
}

int estimateMaxNGhostPerCornerNeighbor(double rho, double cellsize[DIM])
{
 /* Function to estimate how many particle this region will get from a
    neighbor as ghost particles. Based on the size of a cell, and the
    average system density.  The formula uses a base of 2 particles
    (for no reason), then allows for 55% systematic variations in the
    density and for sqrt(N) statistical fluctuations in the number of
    particles. */
    double nonfluctuatingN = 2.0 + 1.55*cellsize[0]*cellsize[1]*cellsize[2]*rho;
    return (int)(nonfluctuatingN + 2.0*sqrt(nonfluctuatingN));
}

int estimateNmax(int Nmin, double rho, double localL[DIM], double cellsize[DIM])
{
 /* Estimate estimate of the maximum number of particles a given
    processor will have to store at most (including ghost
    particles) */
    int result = estimateMaxNSelf(Nmin, rho, localL) + 6*estimateMaxNGhostPerFaceNeighbor(rho, localL, cellsize)
        + 12*estimateMaxNGhostPerEdgeNeighbor(rho, localL, cellsize) + 8*estimateMaxNGhostPerCornerNeighbor(rho, cellsize);
    return result;
}

long long estimateNpairsmax(int Nmin, double rho, double localL[DIM], double cellsize[DIM])
{
 /* Estimate estimate of the maximum number of interactions a given
    processor will have to compute per time step (including with ghost
    particles) */
    /* based on the density estimate the following */
    /* 1. the maximum number of particle in the region*/
    long long Nmaxself = estimateMaxNSelf(Nmin, rho, localL);
    /* 2. the maximum number of neighbor-region particles that will
       interact with this region.*/
    long long ImaxPerParticle = (int)(3 + 26*cellsize[0]*cellsize[1]*cellsize[2]*rho);
    /* note: these should be fairly conservative estimates, but they
       are still estimate. */
    return Nmaxself*ImaxPerParticle;
}

