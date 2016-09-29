#ifndef ESTIMATESH
#define ESTIMATESH

/** 
 * Function to estimate how many particle this region will hold at
 * most. Based on the volume of the region and its density. Must be
 * larger than Nmin.  The formula uses a base of 13 particles (for no
 * reason), then allows for 2% systematic variations in the density
 * and for 2*sqrt(N) statistical fluctuations in the number of
 * particles. 
 *
 * \param Nmin minimum number of particles in region
 * \param rho density of the system
 * \param localL array of DIM doubles with the size of this region
 *
 * \result Estimate of how many particle this region will hold at
 *         most.
 */
int estimateMaxNSelf(int Nmin, double rho, double localL[DIM]);

/** Function to estimate how many particle this region will get from
 * an edge neighbor as ghost particles. Based on the maximal surface
 * of an edge of the region, the size of a cell, and the average
 * system density.  The formula uses a base of 6 particles (for no
 * reason), then allows for 55% systematic variations in the density
 * and for sqrt(N) statistical fluctuations in the number of
 * particles.
 *
 * \param rho density of the system
 * \param localL size of the system in three dimensions
 * \param cellsize size of a cell in three dimensions
 *
 * \result estimate how many particle this region will get from a face neighbor.
 */
int estimateMaxNGhostPerFaceNeighbor(double rho, double localL[DIM], double cellsize[DIM]);

/**
 * Function to estimate how many particle this region will get from an
 * face neighbor as ghost particles. Based on the temperature, the
 * time step, and the range of the potential, as well as the maximal
 * surface of an edge of the region.  There is almost no basis to this
 * formulat. It uses a base of 1200 particles, scales with temperature
 * and time step as a square root, allows with sqrt dependency of the
 * result.  a*a is used as an estimate of the collisional cross
 * section.  This is inherently the shakiest estimate of the bunch, as
 * initial conditions can lead to a large number of exchanged
 * particles initially.
 */
int estimateMaxNCrossThroughFace(double rho, double localL[DIM], double T, double dt, double a);

/* 
 * Function to estimate how many particle this region will get from a
 * neighbor as ghost particles. Based on the maximal length of an edge
 * of the region, the size of a cell, and the average system density.
 * The formula uses a base of 6 particles (for no reason), then allows
 * for 55% systematic variations in the density and for sqrt(N)
 * statistical fluctuations in the number of particles. 
 *
 * \param rho density of the system
 * \param localL size of the system in three dimensions
 * \param cellsize size of a cell in three dimensions
 *
 * \result estimate how many particle this region will get as ghost particles.
 */
int estimateMaxNGhostPerEdgeNeighbor(double rho, double localL[DIM], double cellsize[DIM]);

/**
 * Function to estimate how many particle this region will get from a
 * neighbor as ghost particles. Based on the size of a cell, and the
 * average system density.  The formula uses a base of 2 particles
 * (for no reason), then allows for 55% systematic variations in the
 * density and for sqrt(N) statistical fluctuations in the number of
 * particles. 
  *
  * \param rho density of the system
  * \param cellsize size of a cell in three dimensions
  *
  * \result estimate how many particle this region will get from a face neighbor.

 */
int estimateMaxNGhostPerCornerNeighbor(double rho, double cellsize[DIM]);

/**
 * Estimate estimate of the maximum number of particles a given
 * processor will have to store at most (including ghost particles)
  *
 * \param Nmin minimum number
  * \param rho density of the system
  * \param localL size of the system in three dimensions
  * \param cellsize size of a cell in three dimensions
  *
  * \result estimate how many particle this region will get from a face neighbor.

 */
int estimateNmax(int Nmin, double rho, double localL[DIM], double cellsize[DIM]);

/**
 * Estimate estimate of the maximum number of interactions a given
 * processor will have to compute per time step (including with ghost
 * particles)
 *
 * \param Nmin minimum number
 * \param rho density of the system
 * \param localL size of the system in three dimensions
 * \param cellsize size of a cell in three dimensions
 *
 * \result estimate of the maximum number of interactions a given
 *         processor will have to compute per time step.
 */
long long estimateNpairsmax(int Nmin, double rho, double localL[DIM], double cellsize[DIM]);

#endif
