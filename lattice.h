#ifndef LATTICEH
#define LATTICEH

#define NO_LATTICE_POSITION -1

typedef struct lattice_t {
    double    a;         /**<lattice distance (cubic) */
    double    L;         /**<length of system in each direction */
    int       i;         /**<indices of last point, in each direction */
    int       j;
    int       k;
    long long maxpoints; /**<maximum number of points */
    long long maxperdim; /**<maximum number of points per direction */
    int       starti;    /**<first point in i direction */
    int       startj;    /**<first point in j direction */
    int       startk;    /**<first point in k direction */
    int       stopi;     /**<last point in i direction */
    int       stopj;     /**<last point in j direction */
    int       stopk;     /**<last point in k direction */
} lattice_t;

/** Initialize a lattice point generator for a system of size LxLxL
   with N points. Coordinates of the points will lie between -L/2 and
   +L/2 in each direction.*/
lattice_t init_lattice(double L, long long N);

/** Initialize a lattice point generator for a part of a system of size
   LxLxL with total number of N points (with points between -L/2 and
   L/2). Part is defined by its origin and dimensions
   localL. Coordinates of the points will lie between origin[d] and
   origin[d]n+localL[d] in each direction d=0,1,2. */
lattice_t init_partial_lattice(double L, long long N, double* origin, double* localL);

/** Fill x, y, z with the next point of the lattice. */
long long make_lattice_position(lattice_t* lat, double* x, double* y, double* z);

#endif
