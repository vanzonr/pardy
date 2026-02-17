///
/// @file cells.cpp
///
/// @brief TDB
///
/// @author Ramses van Zon
/// @date 2026
///
#include "cells.h"
#include "global.h"
#include "debug.h"

static bool order_by_atom_super_cell(const atom_t& p1, const atom_t& p2)
{
    return p1.c < p2.c;
}

void cellDivide(rvector<atom_t>& atoms, system_t& sys)
{
    /* Setup cell structure */
    /* There is a global (or 'tot') number of cells, which is set
       by the total system size L and the true cell size. The true
       cell size is larger or equal to the minimum cell size that
       is given as the min_cell_size argument to this function,
       such that the system size is a multiple of the true cell
       size.

       What is actually used in this cell division is the 'local'
       number of cells, which is derived from inspecting which
       cells contain particles, and using the smallest (3d)
       rectangular section of cells. This reduces the number of
       cells to consider (or at most keeps it the same), but will
       become important in the mpi version, where each process
       will contain only a subset of the particles.
    */
    #pragma omp parallel for 
    for (int i = 0; i < sys.N; ++i) {
        int cx = (atoms[i].rx + 0.5*sys.L)/sys.cellsize[0];
        int cy = (atoms[i].ry + 0.5*sys.L)/sys.cellsize[1];
        int cz = (atoms[i].rz + 0.5*sys.L)/sys.cellsize[2];
        /* check that all atoms are where they should be: */
        assert_le(sys.minc[0], cx); assert_le(cx, sys.maxc[0]);
        assert_le(sys.minc[1], cy); assert_le(cy, sys.maxc[1]);
        assert_le(sys.minc[2], cz); assert_le(cz, sys.maxc[2]);
        /* if run without debugging, at least do this: */
        if (cx < sys.minc[0])
            cx = sys.minc[0];
        if (cx > sys.maxc[0])
            cx = sys.maxc[0];
        if (cy < sys.minc[1])
            cy = sys.minc[1];
        if (cy > sys.maxc[1])
            cy = sys.maxc[1];
        if (cz < sys.minc[2])
            cz = sys.minc[2];
        if (cz > sys.maxc[2])
            cz = sys.maxc[2];
        atoms[i].cx = cx;
        atoms[i].cy = cy;
        atoms[i].cz = cz;
    }
    /* the cell list is an ordered array of particles indices such
       that the indices of particles in the same cell are together. We
       need to know where in that list a particular cell s starts, and
       how many particles are in that cell. */
    #pragma omp parallel for
    for (int i = 0; i < sys.N; ++i)
        atoms[i].c = compute_super_cell_index(atoms[i].cx, atoms[i].cy, atoms[i].cz, sys.minc, sys.nc);   
    for (int i = 0; i < sys.ncprod; ++i)
        sys.n_in_cell[i] = 0;
    for (int i = 0; i < sys.N; ++i) 
        sys.n_in_cell[atoms[i].c] += 1;
    /* sort so that the cells form */
    std::sort(atoms.begin(), atoms.begin() + sys.N, order_by_atom_super_cell);
    /* find where each cell starts within the atom list */
    size_t current = 0;
    for (int i = 0; i < sys.ncprod; ++i) {
        sys.start_of_cell[i] = current;
        current += sys.n_in_cell[i];
    }
}

void findPairsFromOtherCells(interaction_pairs_t& p, system_t& sys,
                             int ai, int cax, int cay, int caz, int nneighbor_cells)
{
    #define max_nneighbor_cells 26
    const vec<tiny_int> neighbor_cell_offset[max_nneighbor_cells]
        = {   { 0,  1, 0}, { 1, -1, 0}, { 1, 0, 0}, {1, 1, 0},
              {-1, -1, 1}, {-1,  0, 1}, {-1, 1, 1},
              { 0, -1, 1}, { 0,  0, 1}, { 0, 1, 1},
              { 1, -1, 1}, { 1,  0, 1}, { 1, 1, 1},
              { 0, -1, 0}, {-1,  1, 0}, {-1, 0, 0}, {-1,-1, 0},
              { 1,  1,-1}, { 1,  0,-1}, { 1,-1,-1},
              { 0,  1,-1}, { 0,  0,-1}, { 0,-1,-1},
              {-1,  1,-1}, {-1,  0,-1}, {-1,-1,-1} };
    size_t k = p.npairs;
    for (int dc = 0; dc < nneighbor_cells; ++dc) {
        /* cbx is the integer x coordinate of cell b. db contains the
           integer offsets to apply to particles in that cell (times L
           for actual offset) */
        int cbx = cax + neighbor_cell_offset[dc][0];
        packed_int db = 0;
        if  (cbx < 0) {
            /* we have a negative cell; this needs to be
               mapped to the cell true number as stored by
               adding totncx. For particles in that cell, they
               are stored in that true cell as the coordinates
               they have there; we will have to put them back
               a distance -L to count as neighbors. */
            cbx += sys.totnc[0]; 
            db |= minux; /* puts them back -L */
        } else if (cbx >= sys.totnc[0]) {
            cbx -= sys.totnc[0];
            db |= plusx;
        } 
        if (cbx < sys.minc[0] || cbx > sys.maxc[0])
            continue; /* skip the cell if we don't have those particle; someone else's reponsibility. */
        /* y and z directiion same way */
        int cby = cay + neighbor_cell_offset[dc][1];
        if (cby < 0) {
            cby += sys.totnc[1];
            db |= minuy;
        } else if (cby >= sys.totnc[1]) {
            cby -= sys.totnc[1];
            db |= plusy;
        }
        if (cby < sys.minc[1] || cby > sys.maxc[1])
            continue; /* skip the cell if we don't have those particle; someone else's reponsibility. */
        int cbz = caz + neighbor_cell_offset[dc][2];
        if (cbz < 0) {
            cbz += sys.totnc[2];
            db |= minuz;
        } else if (cbz >= sys.totnc[2]) {
            cbz -= sys.totnc[2];
            db |= plusz;
        }
        if (cbz < sys.minc[2] || cbz > sys.maxc[2])
            continue; /* skip the cell if we don't have those particle; someone else's reponsibility. */
        /* now consider all particles in cell b */
        const int cb = compute_super_cell_index(cbx, cby, cbz, sys.minc, sys.nc);
        const int bjbegin = sys.start_of_cell[cb];
        const int bjend = sys.start_of_cell[cb] + sys.n_in_cell[cb];
        for (int bj = bjbegin; bj < bjend; ++bj) {
            p.pairi[k-bjbegin+bj] = ai;
            p.pairj[k-bjbegin+bj] = bj;
            p.dj[k-bjbegin+bj] = db;            
        }
        k += bjend-bjbegin;
    }
    p.npairs = k;
}

void findGhostPairsFromCells(interaction_pairs_t& p, int ghost_count, rvector<atom_t>& ghost_atoms, int ghost_offset, system_t& sys)
{
    /* Determines the interaction pairs between system particles
     * stored in cells with ghost_count ghost_atoms.  IN: cells,
     * ghost_count, ghost_atoms, ghost_offset OUT: npairs, pairi,
     * pairj, djx, djy, djz . The values of pairi are ghost particles,
     * shifted by ghost_offset, while pairj are system particles */
    /* Relative neighbor cell coordinates. */
    p.npairs = 0; 
    for (int ai = 0; ai < ghost_count; ++ai) {
        int cax = (ghost_atoms[ai].rx + 0.5*sys.L)/sys.cellsize[0];
        int cay = (ghost_atoms[ai].ry + 0.5*sys.L)/sys.cellsize[1];
        int caz = (ghost_atoms[ai].rz + 0.5*sys.L)/sys.cellsize[2];
        findPairsFromOtherCells(p, sys, ai + ghost_offset, cax, cay, caz, 26);
    }
}

void findPairsFromCells(interaction_pairs_t& p, system_t& sys)
{
    /* Collect the pairs of particles in the same or neighboring cells. 
       Output are p.npairs, p.pairi[], p.pairj[] and
       p.dj{x,y,z}[], where i=p.pairi[k] and j=p.pairj[k] are
       neighbors (k=0..*npairs-1), and the integer values *dj{x,y,z}
       indicate if the neighborship is realized through the boundary
       conditions, such that the neighbor position for particle j
       should be taken as atom[j].rx + djx*Lx (and similar for y and
       z). */

    /* Relative neighbor cell coordinates. Note that while there are
       26 neighbor cells, we only need to consider 13 as the others
       are reached from the other cell. */
    p.npairs = 0;
    for (int ca = 0; ca < sys.ncprod; ++ca) {
        int cax = super_cell_index_2_cx(ca, sys.minc, sys.nc);
        int cay = super_cell_index_2_cy(ca, sys.minc, sys.nc);
        int caz = super_cell_index_2_cz(ca, sys.minc, sys.nc);
        const int aibegin = sys.start_of_cell[ca];
        const int aiend = sys.start_of_cell[ca] + sys.n_in_cell[ca];
        for (int ai = aibegin; ai < aiend; ++ai) {
            /* now the task is to find which particles it interacts with 
                first, within the same cell (note the start index!)*/
            long long k = p.npairs;
            for (int aj = ai+1; aj < aiend; ++aj) {
                p.pairi[k] = ai;
                p.pairj[k] = aj;
                p.dj[k] = 0;
                ++k;
            }
            p.npairs = k;
            /* particles in other cells */
            findPairsFromOtherCells(p, sys, ai, cax, cay, caz, 13);
        }
    }
}
