#include "lattice.h"
#include <cmath>

lattice_t init_lattice(double L, long long N)
{
 /* Initialize a lattice point generator for a system of size LxLxL
    with N points. Coordinates of the points will lie between -L/2 and
    +L/2 in each direction.*/
    const int Lovera = (int)(cbrt(1.0*N)+0.99999999999);
    const double a = L/Lovera; /* lattice distance */
    return (lattice_t){a, L, NO_LATTICE_POSITION, NO_LATTICE_POSITION, NO_LATTICE_POSITION,
            N, Lovera, 0, 0, 0, Lovera-1, Lovera-1, Lovera-1, };
}

lattice_t init_partial_lattice(double L, long long N, const vec<double>& origin, const vec<double>& localL)
{
 /* Initialize a lattice point generator for a part of a system of size
    LxLxL with total number of N points (with points between -L/2 and
    L/2). Part is defined by its origin and dimensions
    localL. Coordinates of the points will lie between origin[d] and
    origin[d]n+localL[d] in each direction d=0,1,2. */
    const double releps = 1.0e-5; /* needed for floating point tolerance */
    const int Lovera = (int)(cbrt(1.0*N)+0.99999999999);
    const double a = L/Lovera; /* lattice distance */
    int start[3], stop[3];
    for (int d=0;d<3;d++) {
        start[d] = (origin[d]+0.5*L)/a;
        while (start[d]*a < origin[d]+0.5*L)
            start[d]++;
        stop[d] = (origin[d]+0.5*L+localL[d])/a+1;
        while (stop[d]*a*(1+releps) >= origin[d]+0.5*L+localL[d] ||
               stop[d]*a*(1-releps) >= origin[d]+0.5*L+localL[d]   ) {
            stop[d]--;
        }
    }
    return (lattice_t){a, L, NO_LATTICE_POSITION, NO_LATTICE_POSITION, NO_LATTICE_POSITION,
            N, Lovera, start[0], start[1], start[2], stop[0], stop[1], stop[2]};
}

long long make_lattice_position(lattice_t& lat, double& x, double& y, double& z)
{
    /* Given a lattice_t initialized with init_lattice_*, fills x, y, z
       with the next particle on the lattice. It returns the index
       of that particle, or -1 if there are no more particles.
    */
    if (lat.i == NO_LATTICE_POSITION ||
        lat.j == NO_LATTICE_POSITION ||
        lat.k == NO_LATTICE_POSITION)
    {
        lat.i = lat.starti;
        lat.j = lat.startj;
        lat.k = lat.startk;
    } else {
        ++(lat.i);
        if (lat.i > lat.stopi) {              /* outside box? */
            lat.i = lat.starti;               /* then put back */
            ++(lat.j);                         /* and move to next y level */
            if (lat.j > lat.stopj) {          /* outside box? */
                lat.j = lat.startj;           /* then put back */
                ++(lat.k);                     /* and move to next z level */
                if (lat.k > lat.stopk)        /* outside box? */
                    return NO_LATTICE_POSITION; /* no more points */
            }
        }
    }
    long long index = (lat.k*lat.maxperdim+lat.j)*lat.maxperdim+lat.i;
    if (index >= lat.maxpoints)
        return NO_LATTICE_POSITION; /* no more points */
    else {
        x = lat.i*lat.a - 0.5*lat.L;
        y = lat.j*lat.a - 0.5*lat.L;
        z = lat.k*lat.a - 0.5*lat.L;
        return index;
    }
}
