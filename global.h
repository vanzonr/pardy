/* a few global setting for all of pardy */
#ifndef GLOBALH
#define GLOBALH
/* global_rank and global_size should be set, once, at the start of main. */
extern int global_size; 
extern int global_rank;
const int global_root = 0;
#define i_am_root (global_rank == global_root)
#define DIM 3
typedef signed char tiny_int;

#endif
