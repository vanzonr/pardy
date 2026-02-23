#ifndef _LCG2_H_
#define _LCG2_H_

#include <memory>

class rand48lcg;
class predictable_normal_distribution;
using lcg_t = std::shared_ptr<rand48lcg>;
using lcg_dist_t = std::shared_ptr<predictable_normal_distribution>;

lcg_t lcg_init_rand48(long seed);

lcg_dist_t lcg_init_normal(double mean, double stddev, bool precise=false);

double lcg_dist_apply(lcg_t& lcg, lcg_dist_t& dist);

void lcg_skip(lcg_t& lcg, long n);

#endif
