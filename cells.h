/* cells.h */
#ifndef CELLSH
#define CELLSH

#include "atom.h"
#include "system.h"
#include "forces.h"
#include <rarray>

/**
 * Function to compute a super-cellindex out of the cell indices cx,
 * cy and cz. This super cell index is used to sort the particles into
 * cells, and as an ordering for cells. cx is assumed to lie between
 * minc[0] and minc[0]+nc[0], cy between minc[1] and minc[1]+nc[1],
 * and cz between minc[2] and minc[2]+nc[2].  The equation is:
 * super-index = ((z-minc[2])*nc[1]+(y-minc[1]))*nc[0]+(x-minc[0]);
 *
 * \param cx cell index in the x direction.
 * \param cy cell index in the y direction.
 * \param cz cell index in the z direction.
 * \param minc[] is an array of 3 elements which are the minimium cell indices in each direction.
 * \param nc[] is an array of 3 elements which are the number of cells in each direction.
 *
 * \result the super index
 */
static inline int compute_super_cell_index(int x, int y, int z, const vec<int>& minc, const vec<int>& nc)
{
    return ((z-minc[2])*nc[1]+(y-minc[1]))*nc[0]+(x-minc[0]);
}

/**
 * Function to extract the x cell index.
 *
 * \param c the super index
 * \param minc[] is an array of 3 elements which are the minimium cell indices in each direction.
 * \param nc[] is an array of 3 elements which are the number of cells in each direction.
 *
 * \result the cell index in the x directions.
 *
 * \seealso compute_super_cell_index
 */
static inline int super_cell_index_2_cx(int c, const vec<int>& minc, const vec<int>& nc)
{
    return (c%(nc[0]))+minc[0];
}

/**
 * Function to extract the y cell index.
 *
 * \param c the super index
 * \param minc[] is an array of 3 elements which are the minimium cell indices in each direction.
 * \param nc[] is an array of 3 elements which are the number of cells in each direction.
 *
 * \result the cell index in the y directions.
 *
 * \see compute_super_cell_index
 */
static inline int super_cell_index_2_cy(int c, const vec<int>& minc, const vec<int>& nc)
{
    return ((c/nc[0])%(nc[1]))+minc[1];
}

/**
 * Function to extract the z cell index.
 *
 * \param c the super index
 * \param minc[] is an array of 3 elements which are the minimium cell indices in each direction.
 * \param nc[] is an array of 3 elements which are the number of cells in each direction.
 *
 * \result the cell index in the z directions.
 *
 * \see compute_super_cell_index
 */
static inline int super_cell_index_2_cz(int c, const vec<int>& minc, const vec<int>& nc)
{
    return (c/(nc[0]*nc[1]))+minc[2];
}

/**
 * Setup cell structure 
 *
 * \param atoms the atoms
 * \param the system to subdivide in cells and organize particles.
 */
void cellDivide(rvector<atom_t>& atoms, system_t& sys);

/**
 * Function to find pairs with particles in neighboring cells given a
 * cell.
 *
 * \param p a reference to an interaction_pairs_t structure that will
 *          hold the interaction pairs, i.e., in the indices of the
 *          pair of particles that interact.
 * \param sys a reference to the system 
 * \param ai  the index of the first particle
 * \param cax the x-cell index of the cell in which particle ai resides
 * \param cay the x-cell index of the cell in which particle ai resides
 * \param caz the x-cell index of the cell in which particle ai resides
 * \param nneighbor_cells the number of neighbor cells to consider
 */
void findPairsFromOtherCells(interaction_pairs_t& p, system_t& sys,
                             int ai, int cax, int cay, int caz, int nneighbor_cells);

/**
 * Function to find pairs with particles in neighboring cells given a
 * set of ghost particles.
 *
 * \param p a reference to an interaction_pairs_t structure that will
 *          hold the interaction pairs, i.e., in the indices of the
 *          pair of particles that interact.
 * \param ghost_count the number of ghost particles
 * \param ghost_atoms the ghost particles
 * \param ghost_offset a number to add to the index of the ghost atoms
 *         in the ghost_atoms arrays when constructing the pair list.
 * \param sys a reference to the system 
 */
void findGhostPairsFromCells(interaction_pairs_t& p, int ghost_count, rvector<atom_t>& ghost_atoms, int ghost_offset, system_t& sys);

/** 
 * The function 'findPairsFromCells collects the pairs of particles
 * in the same or neighboring cells.  Output are p->npairs,
 * p->pairi[], p->pairj[] and p->dj{x,y,z}[], where i=p->pairi[k] and
 * j=p->pairj[k] are neighbors (k=0..*npairs-1), and the integer
 * values p->dj{x,y,z}[k] indicate if the neighborship is realized through
 * the boundary conditions, such that the neighbor position for
 * particle j should be taken as atom[j].rx + p->djx[j]*Lx (and similar for
 * y and z). 
 *
 * \param p a reference to an interaction_pairs_t structure that will
 *          hold the interaction pairs, i.e., in the indices of the
 *          pair of particles that interact.
 * \param sys a reference to the system 
 */
void findPairsFromCells(interaction_pairs_t& p, system_t& sys);

#endif
