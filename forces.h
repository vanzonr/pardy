#ifndef FORCESH
#define FORCESH

#include "global.h"
#include "atom.h"
#include "parallelwork.h"
#include <rarray>

static const double rc  = 3.0;   /**< outer cutoff radius */
static const double rcp = 2.5;   /**< inner cutoff, or start of smoothing */

typedef signed char packed_int;
#define plusx 0x01
#define minux 0x02
#define plusy 0x04
#define minuy 0x08
#define plusz 0x10
#define minuz 0x20
#define pack2x(pack)  (((pack)&minux)?-1:(((pack)&plusx)?1:0))
#define pack2y(pack)  (((pack)&minuy)?-1:(((pack)&plusy)?1:0))
#define pack2z(pack)  (((pack)&minuz)?-1:(((pack)&plusz)?1:0))

struct interaction_pairs_t {
    long long   npairs; /**< number of currently stored interaction pairs */
    rvector<int> pairi;  /**< current index in atoms of the first atom of each pair */
    rvector<int> pairj;  /**< current index in atoms of the second atom of each pair */
    rvector<packed_int> dj;     /**< shifts in x,y, and z to apply to second
                             atom (either -1, 0, or 1); these shifts
                             accounts for periodic boundary conditions */
};

/**
 *
 */
void interaction_pairs_alloc(interaction_pairs_t& p, long long Npairsmax);

/**
 *
 */
void interaction_pairs_free(interaction_pairs_t& p);

/**
 *
 */
double computeForces(int N, rvector<atom_t>& atoms, interaction_pairs_t& p, double L, bool accum, parallel_work_t& work);

#endif
