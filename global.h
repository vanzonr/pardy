/* a few global setting for all of pardy */
#ifndef GLOBALH
#define GLOBALH
/* global_rank and global_size should be set, once, at the start of main. */
extern int global_size; 
extern int global_rank;
#define global_root 0
#define i_am_root (global_rank == 0)
#define DIM 3
typedef signed char tiny_int;

/**
 * Helper function to find maximum of three double precision numbers 
*/
static inline double max3(double x1, double x2, double x3)
{
    double m12 = ((x1>x2) ? x1 : x2);
    return ((m12>x3) ? m12 : x3);
}
#endif
