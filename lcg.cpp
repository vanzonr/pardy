#include "lcg48.h"
#include "lcg.h"

lcg_t lcg_init_rand48(long seed)
{
    return lcg_t(new rand48lcg(seed));
}

lcg_dist_t lcg_init_normal(double mean, double stddev, bool precise)
{
    return lcg_dist_t(new predictable_normal_distribution(mean, stddev, precise));
}

double lcg_dist_apply(lcg_t& lcg, lcg_dist_t& dist)
{
    return (*dist)(*lcg);
}

void lcg_skip(lcg_t& lcg, long n)
{
    (*lcg).discard(n);
}
