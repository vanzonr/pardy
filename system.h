///
/// @file system.h
///
/// @brief Defines  a struct for system parameters and its MPI type.
///
/// @author Ramses van Zon
/// @date 2026
///
#ifndef SYSTEMH
#define SYSTEMH
#include <mpi.h>
#include <cstddef>
#include "global.h"
#include <rarray>

struct system_t {
    long long Ntot;      /**< number of particles in the full domain (parameter) */
    double    rho;       /**< density of the system (parameter)*/
    double    T0;        /**< initial temperature (parameter) */
    double    runtime;   /**< how long to run (parameter) */
    double    dt;        /**< duration of single time step (parameter) */
    long      seed;      /**< global random seed (parameter) */
    double    equil;     /**< when to start measuring energies etc. (parameter)*/
    bool      usecells;  /**< 1 if cells are to be used (parameter) */
    double    L;         /**< length of the box (derived parameter) */  
    int       N;         /**< number of particles in the subdomain*/
    double    K;         /**< kinetic energy in the subdomain*/
    double    U;         /**< potential energy in the subdomain */
    /**< for mpi domain decomposition: */
    MPI_Comm  comm;      /**< communicator (cartesian variant of COMM_WORLD) */
    int       nprocs;    /**< total number of mpi processes */
    vec<int>  np;        /**< size of the 3d processes grid */
    int       rank;      /**< rank of this process in communicator comm */
    int       neighborrank[27]; /**< ranks of the neighbors. Entry 13 is
                            this process, 13+/-1 is the process to the
                            forward or backward in the x direction,
                            13+/-3 is the process right or left in the
                            y direction, 13+/-9 is the process up or
                            down in the z direction. */ 
    bool      neighborsendblock[27]; /**< which of the above ranks to
                            skip when sending ghost particles, because
                            it would deliver the same ghost data to the
                            same process twice (or more). This can
                            happen if the number of processors in a
                            given dimension is two, and that direction
                            only allows for 1 cell per
                            processor. Likewise, if the number of
                            processors in a given dimension is one, no
                            communication is needed in that
                            direction. */
    vec<double> localL;  /**< size of the local subdomain for this mpi process */
    vec<double> origin;  /**< origin of the local subdomain for this mpi process */
    /**< for cell structure: */
    vec<int>    totnc;   /**< total number of cells in the whole system */
    vec<int>    nc;      /**< number of cells in the subdomain in each direction */
    int         ncprod;  /**< product of nc[0..2], i.e., total number of cells in this subdomain */
    vec<int>    minc;    /**< first cell index in 3 dimensions in the subdomain */
    vec<int>    maxc;    /**< last cell index in 3 dimensions in the subdomain */
    vec<double> cellsize;/**< size of the cells in each direction */
    rarray<int,1> start_of_cell; /**< position in atom[] of the first particles in each cell (1d using a super index) */
    rarray<int,1> n_in_cell; /**< number of particles in each cell (1d using a super index) (i.e., subsequent n_in_cell-1 particles are also in that cell */
};

extern MPI_Datatype MPI_PARAMETERS;

bool sanity_check(system_t& sys, double rc, bool i_am_root);
MPI_Datatype define_MPI_PARAMETERS();

/* what is where in neighborranks ? */
#define TOTHEFRONT  (+1) /* x incr */
#define TOTHEBACK   (-1) /* x decr */
#define TOTHERIGHT  (+3) /* y incr */
#define TOTHELEFT   (-3) /* y decr */
#define TOTHETOP    (+9) /* z decr */
#define TOTHEBOTTOM (-9) /* z incr */
#define CENTER 13
#define FRONT       (CENTER+TOTHEFRONT)
#define BACK        (CENTER+TOTHEBACK)
#define RIGHT       (CENTER+TOTHERIGHT)
#define LEFT        (CENTER+TOTHELEFT)
#define TOP         (CENTER+TOTHETOP)
#define BOTTOM      (CENTER+TOTHEBOTTOM)
#define FRONTRIGHT  (CENTER+TOTHEFRONT+TOTHERIGHT)
#define FRONTLEFT   (CENTER+TOTHEFRONT+TOTHELEFT)
#define BACKRIGHT   (CENTER+TOTHEBACK+TOTHERIGHT)
#define BACKLEFT    (CENTER+TOTHEBACK+TOTHELEFT)
#define FRONTTOP    (CENTER+TOTHEFRONT+TOTHETOP)
#define FRONTBOTTOM (CENTER+TOTHEFRONT+TOTHEBOTTOM)
#define BACKTOP     (CENTER+TOTHEBACK+TOTHETOP)
#define BACKBOTTOM  (CENTER+TOTHEBACK+TOTHEBOTTOM)
#define RIGHTTOP    (CENTER+TOTHERIGHT+TOTHETOP)
#define RIGHTBOTTOM (CENTER+TOTHERIGHT+TOTHEBOTTOM)
#define LEFTTOP     (CENTER+TOTHELEFT+TOTHETOP)
#define LEFTBOTTOM  (CENTER+TOTHELEFT+TOTHEBOTTOM)
#define FRONTTOPRIGHT    (CENTER+TOTHEFRONT+TOTHETOP+TOTHERIGHT)
#define FRONTTOPLEFT     (CENTER+TOTHEFRONT+TOTHETOP+TOTHELEFT)
#define FRONTBOTTOMRIGHT (CENTER+TOTHEFRONT+TOTHEBOTTOM+TOTHERIGHT)
#define FRONTBOTTOMLEFT  (CENTER+TOTHEFRONT+TOTHEBOTTOM+TOTHELEFT)
#define BACKTOPRIGHT     (CENTER+TOTHEBACK+TOTHETOP+TOTHERIGHT)
#define BACKTOPLEFT      (CENTER+TOTHEBACK+TOTHETOP+TOTHELEFT)
#define BACKBOTTOMRIGHT  (CENTER+TOTHEBACK+TOTHEBOTTOM+TOTHERIGHT)
#define BACKBOTTOMLEFT   (CENTER+TOTHEBACK+TOTHEBOTTOM+TOTHELEFT)

#endif
