///
/// @file lcg.h
///
/// @brief Implementation of a linear congruential generator with
/// log(n) skipping. Able to mimic the rand48 family of random number
/// functions from C.
///
/// USAGE
/// 
/// The linear congruential generator requires three parameters (a, b,
/// and m), and a state variable (s), which are kept in a 'lcg_t'
/// structure.
///
/// FUNCTIONS
///
///  lcg_t lcg_init(uint64_t a, uint64_t b, uint64_t m, uint64_t seed);
///  lcg_t lcg_init_rand48(long int seed);
///
///   The function lcg_init() can be used to set these variable and the
///   'seed' for the state.  Default values for a, b and m to mimic the
///   behaviour of the 'rand48' family functions are provided through
///   the function lcg_init_rand48(), which only requires a long int
///   'seed' parameter, which in the same way as srand48 would
///   have. Both lcg_init functions return a 'lcg_t' to be used in
///   calls to the other lcg_* functions.
///
///  void lcg_seed(lcg_t* lcg, uint64_t seed);
///
///   The function lcg_seed() resets the seed of the random number
///   generator 'lcg'.
///
///  float     lcg_float (lcg_t* lcg);
///  double    lcg_double(lcg_t* lcg);
///  int32_t   lcg_int32 (lcg_t* lcg);
///  int32_t   lcg_pint32(lcg_t* lcg);
///  uint32_t  lcg_uint32(lcg_t* lcg);
///  uint64_t  lcg_int64(lcg_t* lcg);
/// 
///   These functions return the next random number for the random
///   number generator 'lcg', in the form of a float, double, 32-bit
///   integer, and 64-bit integer, respectively. The three 32-bit
///   versions differ in the range of returned values, according to the following table:
///   -----------  ---------------
///   function     range
///   -----------  ---------------
///   lcg_int32    [-2^31,2^31-1]
///   lcg_pint32   [0,2^31-1]
///   lcg_uint32   [0,2^32-1]
///   -----------  ---------------
///
///  void lcg_skip(lcg_t* lcg, int64_t n);
/// 
///   The function lcg_skip() skips the next 'n' random numbers, in
///   O(log n) time. In the current implementation, if n is negative,
///   the function will abort (or do nothing if compiled with
///   NDEBUG). The skip-ahead feature can be useful in parallel
///   processing.
/// 
///  void   lcg_srand48(lcg_t* lcg, long int seed);
///  long   lcg_lrand48(lcg_t* lcg);
///  long   lcg_mrand48(lcg_t* lcg);
///  double lcg_drand48(lcg_t* lcg);
/// 
///   These functions mimic the behavior of rand48. The 32-bit seed
///   passed to 'lcg_srand48' is treated in the same way to generate
///   the 48-bit seed, and the same algorithms are used to generate
///   long integers ('lcg_lrand48()'), signed long integers
///   ('lcg_mrand48()'), and double precision floating point numbers
///   ('lcg_drand48()').  If 'lcg' is initialized with
///   lcg_init_rand48, and the seed is subsequently set with
///   lcg_srand48, the set of random numbers generated should be
///   identical to using srand48(), lrand48(), mrand48(), and
///   drand48().  The added benefit of using lcg over rand48 is that
///   skip-ahead on 'lcg' is still possible (also rand48 is not in c99
///   or c11 standard).
///
/// IMPLEMENTATION
///
/// Can be found in lcg.cpp.
///
/// @author Ramses van Zon
/// @date 2026
///
#ifndef LCG_HEADER_
#define LCG_HEADER_

#include <stdint.h>

struct lcg_t {
    const uint64_t a;
    const uint64_t b;
    const uint64_t m;
    uint64_t s;
    double x2;
    int have;
};

lcg_t lcg_init(uint64_t a, uint64_t b, uint64_t m, uint64_t seed);
lcg_t lcg_init_rand48(int32_t seed);
lcg_t lcg_init_alt35(int64_t seed);

void lcg_seed(lcg_t& lcg, uint64_t seed);
void lcg_skip(lcg_t& lcg, int64_t n);

float    lcg_float(lcg_t& lcg);
double   lcg_double(lcg_t& lcg);
int32_t  lcg_int32(lcg_t& lcg);
int32_t  lcg_pint32(lcg_t& lcg);
uint32_t lcg_uint32(lcg_t& lcg);
uint64_t lcg_int64(lcg_t& lcg);

double   lcg_normal(lcg_t& lcg);

void     lcg_srand48(lcg_t& lcg, int32_t seed);
long int lcg_lrand48(lcg_t& lcg);
long int lcg_mrand48(lcg_t& lcg);
double   lcg_drand48(lcg_t& lcg);

#endif
