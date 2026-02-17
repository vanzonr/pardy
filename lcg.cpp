/*
 * lcg.cpp - Implementation of a linear congruential generator with
 *           log(n) skipping. For documentation, see the header file
 *           lcg.h.
 *
 * (C) 2016 Ramses van Zon, SciNet HPC Consortium, Toronto, Canada
 */

#include "lcg.h"
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>

const double pi=3.1415926535897932; /* value M_PI not present in c99/c11 */

static void lcg_advance(lcg_t& lcg)
{
    /* Perform one step in the generator. Not exported, as
     * lcg_skip(lcg,1) does the same */
    lcg.s = (lcg.s * lcg.a + lcg.b) % lcg.m;
}

static uint64_t rand2lcgseed(int32_t seed)
{
    /* To initialize the internal buffer r(n) of drand48, lrand48, and
     * mrand48, the 32 bits of the seed value are copied into the
     * upper 32 bits of its state buffer (a "short int r[n]" with
     * n=3), with the lower 16 bits of r[n] arbitrarily being set to
     * 0x330e.  */
    uint64_t seedll = 0;
    memcpy(&seedll, &seed, sizeof(seed));
    return (seedll<<16) + 0x330e;
}

static int32_t u64_2_s32(uint64_t x, int n)
{
    /* The algorithm that mrand48 uses to extract a signed 32-bit
     * value from the 48-bit unsigned long value is to do a 16 bits
     * right-shift. The extra mask is to ensure there is not overflow
     * for the 32-bit long value.  The top bit becomes the sign, so
     * the result may be negative. Here generalized for n bits
     * shifted, rand48 uses n=16.  */
    assert(n >= 0); 
    return (int32_t) ((x>>n)&0xffffffff);
}

static int32_t u64_2_s32_positive(uint64_t x, int n)
{
    /* The algorithm that lrand48 uses to extract a signed long value
     * that is positive from the 48-bit unsigned long value is to do a
     * 16+1 bits right-shift. The +1 in 16+1 is to clear
     * the sign bit, which is further ensured by the 0x7fffffff
     * mask. Here generalized for n+1 bits shifted, rand48 uses
     * n=16.  */
    assert(n >= 0); 
    return (long) ((x>>(n+1)) & 0x7fffffff);
}

static unsigned long u64_2_u32(uint64_t x, int n)
{
    /* An algorithm similar to lrand48 to extract an unsigned long
     * value from the 48-bit unsigned long value is to do a 16 bits
     * right-shift. The extra mask is to ensure there is not overflow
     * for the 32 bit long value. Here generalized for n bits
     * shifted. */
    assert(n >= 0); 
    return (unsigned long) ((x>>n) & 0xffffffff);
}

static int hibit(uint64_t m)
{
    unsigned r = 0;
    while (m >>= 1) 
        r++;
    return r;
}

lcg_t lcg_init(uint64_t a, uint64_t b, uint64_t m, uint64_t seed)
{
    assert(a > 0);
    assert(b >= 0);
    assert(m > 0);
    assert(0 <= seed && seed < m);
    return (lcg_t){a,b,m,seed};
}

lcg_t lcg_init_alt35(int64_t seed)
{
    return (lcg_t){
        .a = 34359738368LL,
        .b = 1LL,
        .m = 1LL << 35,
        .s = (uint64_t)seed
    };
}

lcg_t lcg_init_rand48(int32_t seed)
{
    uint64_t seed64 = rand2lcgseed(seed);
    return (lcg_t){
        .a = 25214903917LL,
        .b = 11LL,
        .m = 1LL << 48,
        .s = seed64
    };
}

void lcg_seed(lcg_t& lcg, uint64_t seed)
{
    assert(lcg.a > 0);
    assert(lcg.b >= 0);
    assert(lcg.m > 0);
    assert(0 <= seed && seed < lcg.m);
    lcg.s = seed;
}

int32_t lcg_int32(lcg_t& lcg)
{
    lcg_advance(lcg);
    int n = hibit(lcg.m)-32;
    return u64_2_s32(lcg.s, n);
}

uint32_t lcg_uint32(lcg_t& lcg)
{
    lcg_advance(lcg);
    int n = hibit(lcg.m)-32;
    return u64_2_u32(lcg.s, n);
}

int32_t lcg_pint32(lcg_t& lcg)
{
    lcg_advance(lcg);
    int n = hibit(lcg.m)-32;
    return u64_2_s32_positive(lcg.s, n);
}

uint64_t lcg_int64(lcg_t& lcg)
{
    lcg_advance(lcg);
    return lcg.s;
}

float lcg_float(lcg_t& lcg)
{
    lcg_advance(lcg);
    return (float)(lcg.s)/(float)(lcg.m);
}

double lcg_double(lcg_t& lcg)
{
    lcg_advance(lcg);
    return (double)(lcg.s)/(double)(lcg.m);
}

/* Inline macro for LCG multiplication, i.e. LGC1 *= LGC2.
 *
 * Based on the formula:
 *
 *  LCG(a1,b1,m) * LCG(a2,b2,m) = LCG((a1*a2)%m, (b1+a1*b2)%m)
 *
 * Note: A macro rather than a function because we want the result
 * in-place. It uses the usual do/while trick for avoiding macro
 * side-effects.
 */
#define LCGMULRIGHT(a1,b1,a2,b2,m)                 \
    do {                                           \
        const uint64_t tmp = a1;                   \
        a1 = (a1*a2)%m;                            \
        b1 = (b1+tmp*b2)%m;                        \
    } while(0)

void lcg_skip(lcg_t& lcg, int64_t n)
{
    /* O(log(n)) skip ahead function for LCG. The algorithm is a
     * variant of modular exponentiation.
     * (https://en.wikipedia.org/wiki/Modular_nonentiation), but with
     * multiplication replaced by multiplication of LCGs (macro
     * above). If n<0, the function will abort through an assert, or,
     * if compile with NDEBUG defined, will not skip at all. */
    assert(n >= 0);
    uint64_t a = 1;
    uint64_t b = 0;
    uint64_t basea = lcg.a % lcg.m;
    uint64_t baseb = lcg.b % lcg.m;
    while (n > 0) {
        if ((n % 2) == 1)
            LCGMULRIGHT(a,b,basea,baseb,lcg.m);
        n >>= 1;
        LCGMULRIGHT(basea,baseb,basea,baseb,lcg.m);
    }
    lcg.s = (lcg.s * a + b) % lcg.m;
}

/* rand48 replacements */

void lcg_srand48(lcg_t& lcg, int32_t seed)
{
    /* srand48 is used to initialize the internal buffer r(n) of
     * drand48, lrand48, and mrand48 such that the 32 bits of the seed
     * value are copied into the upper 32 bits of r(n), with the lower
     * 16 bits of r(n) arbitrarily being set to 0x330e.  */
    lcg_seed(lcg, rand2lcgseed(seed));
}

long lcg_lrand48(lcg_t& lcg)
{
    /* lrand48 and nrand48 return values of type long in the range [0,
     * 2**31-1]. The high-order (31) bits of r(n+1) are loaded into
     * the lower bits of the returned value, with the topmost (sign)
     * bit set to zero. */
    return u64_2_s32_positive(lcg_int64(lcg),16);
}

long lcg_mrand48(lcg_t& lcg)
{
    /* mrand48 and jrand48 return values of type long in the range
     * [-2**31, 2**31-1]. The high-order (32) bits of r(n+1) are
     * loaded into the returned value. */
    return u64_2_s32(lcg_int64(lcg),16);
}

double lcg_drand48(lcg_t& lcg)
{
    /* drand48 and erand48 return values of type double. The full 48
     * bits of r(n+1) are loaded into the mantissa of the returned
     * value, with the exponent set such that the values produced lie
     * in the interval [0.0, 1.0]. */
    union {
        uint64_t long_long_value;
        uint16_t x[4];
    } bitfill;
    bitfill.long_long_value = lcg_int64(lcg);
    const double z = 1.0 / (1L << 16);
    return (z * (z * (z * bitfill.x[0] + bitfill.x[1]) + bitfill.x[2]));
}

double lcg_normal(lcg_t& lcg)
{
    double fac, y1, y2, x1; 
    if (lcg.have == 1) {  /* already one available ? */
        lcg.have = 0;
        return lcg.x2;
    } else {
        /* generate a pair of random variables */
        y1        = lcg_drand48(lcg);
        y2        = lcg_drand48(lcg);
        fac       = sqrt(-2*log(y1));
        lcg.have = 1;
        x1        = fac*sin(2*pi*y2); /* x1 and x2 are now gaussian */
        lcg.x2   = fac*cos(2*pi*y2); /* so store one */
        return x1;              /* and return the other. */
    }
}
