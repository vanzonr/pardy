#ifndef _LCG48_H_
#define _LCG48_H_

#include <cmath>
#include <random>
#include <sstream>
#include <cstring>
#include <cassert>

static uint64_t _rand48_2_lcgseed(int32_t seed)
{
    // To initialize the internal buffer r(n) of drand48, lrand48, and
    // mrand48, the 32 bits of the seed value are copied into the
    // upper 32 bits of its state buffer (a "short int r[n]" with
    // n=3), with the lower 16 bits of r[n] arbitrarily being set to
    // 0x330e. 
    uint64_t seedll = 0;
    memcpy(&seedll, &seed, sizeof(seed));
    return (seedll<<16) + 0x330e;
}

static int32_t _u64_2_s32_positive(uint64_t x, int n)
{
    // The algorithm that lrand48 uses to extract a signed long value
    // that is positive from the 48-bit unsigned long value is to do a
    // 16+1 bits right-shift. The +1 in 16+1 is to clear
    // the sign bit, which is further ensured by the 0x7fffffff
    // mask. Here generalized for n+1 bits shifted, rand48 uses
    // n=16. 
    assert(n >= 0); 
    return (long) ((x>>(n+1)) & 0x7fffffff);
}

static int32_t _u64_2_s32(uint64_t x, int n)
{
    // The algorithm that mrand48 uses to extract a signed 32-bit
    // value from the 48-bit unsigned long value is to do a 16 bits
    // right-shift. The extra mask is to ensure there is not overflow
    // for the 32-bit long value.  The top bit becomes the sign, so
    // the result may be negative. Here generalized for n bits
    // shifted, rand48 uses n=16.
    assert(n >= 0); 
    return (int32_t) ((x>>n)&0xffffffff);
}

template<class UIntType, UIntType a, UIntType c, UIntType m>
class ff_linear_congruential_engine: public std::linear_congruential_engine<UIntType,a,c,m>
{    
  private:
    UIntType get_state() const {
        UIntType x;
        std::stringstream s;
        s << (*this); // get state as a string
        s >> x;       // convert back to unsigned long
        return x;
    }
    static void _lcgmulright(UIntType& a1, UIntType& b1, const UIntType& a2, const UIntType& b2){
        // Based on the formula:
        // LCG(a1,c1,m) * LCG(a2,c2,m) = LCG((a1*a2)%m, (c1+a1*c2)%m)
        const UIntType tmp = a1;
        a1 = (a1*a2)%m;
        b1 = (b1+tmp*b2)%m;
    }
  public:
    typedef typename std::linear_congruential_engine<UIntType,a,c,m>::result_type result_type;
    ff_linear_congruential_engine()
    : std::linear_congruential_engine<UIntType,a,c,m>(std::linear_congruential_engine<UIntType,a,c,m>::default_seed)
    {}
    explicit ff_linear_congruential_engine(result_type value)
    : std::linear_congruential_engine<UIntType,a,c,m>(value)
    {}        
    void discard(uint64_t n) {
        // Faster forward routine for this multiplicative linear
        // congruent generator, O(log n) instead of O(n)
        assert(n >= 0);
        uint64_t s = this->get_state();
        uint64_t newa = 1;
        uint64_t newc = 0;
        uint64_t basea = a % m;
        uint64_t basec = c % m;
        while (n > 0) {
            if ((n % 2) == 1)
                _lcgmulright(newa,newc,basea,basec);
            n >>= 1;
            if (n > 0)
                _lcgmulright(basea,basec,basea,basec);
        }
        s = (s * newa + newc) % m;
        this->seed(s);        
    }
};

class rand48lcg: public ff_linear_congruential_engine<uint64_t,25214903917LL,11ULL,(1ULL<<48)>
{    
  public:
    rand48lcg(long seed)
    : ff_linear_congruential_engine(_rand48_2_lcgseed(seed))
    {}
};

// Random number distribution mimicking lrand, if combined with rand48lcg.
//
// Always produces a positive long between 0 and 2^31-1.
//
// Assumes that the random generator passed into operator() produces
// an unsigned 64-bit integer.
//
class lrand48dist {
  public:
    struct param_type {};
    using result_type = long;    
    template<class URNG>
    result_type operator()(URNG& g) {
        uint64_t number = g();
        return _u64_2_s32_positive(number, 16);
    }
    template<class URNG>
    result_type operator()(URNG& g, const param_type& parm) {
        return this->operator()(g);
    }
    void reset()
    {}
    param_type param() const {
        return param_type{};
    }
    result_type min() const {
        return 0;
    }
    result_type max() const {
        return std::numeric_limits<int32_t>::max();
    }
};

class mrand48dist {
  public:
    struct param_type {};
    using result_type = long;    
    template<class URNG>
    result_type operator()(URNG& g) {
        uint64_t number = g();
        return _u64_2_s32(number, 16);
    }
    template<class URNG>
    result_type operator()(URNG& g, const param_type& parm) {
        return this->operator()(g);
    }
    void reset()
    {}
    param_type param() const {
        return param_type{};
    }
    result_type min() const {
        return std::numeric_limits<int32_t>::min();
    }
    result_type max() const {
        return std::numeric_limits<int32_t>::max();
    }
};

class drand48dist
{
 public:
    struct param_type {};
    using result_type = double;
    template<class URNG>
    result_type operator()(URNG& g) {
        union {
            uint64_t long_long_value;
            uint16_t x[4];            
        } bitfill;
        bitfill.long_long_value = g();
        const double z = 1.0 / (1L << 16);
        return (z * (z * (z * bitfill.x[0] + bitfill.x[1]) + bitfill.x[2]));
    }
    template<class URNG>
    result_type operator()(URNG& g, const param_type& parm) {
        return this->operator()(g);
    }
    void reset()
    {}
    param_type param() const {
        return param_type{};
    }
    result_type min() const {
        return 0.0;
    }
    result_type max() const {
        return 1.0;
    }
};

// seventh-order fit to cos(2*pi*x), accurate up to 1e-5 (but could accumulate)
inline double fastcostwopi(const double x) 
{
    // Derived by fitting values, derivatives and second derivatives
    // at x=0 and x=0.25 (and -by symmetry- also at 0.5, 0.75 and 1.0)
    constexpr const double pi = 3.1415926535897932;
    constexpr const double alpha = 2*pi;
    constexpr const double beta  =   280. -   96.*pi -   2.*pi*pi;
    constexpr const double gamma = -5376. + 1536.*pi +  64.*pi*pi;
    constexpr const double delta = 30720. - 8192.*pi - 512.*pi*pi;
    const double fraction = x - (int64_t)(x);
    const double xi = (fraction > 0.5)?(fraction - 0.75):(0.25 - fraction);
    const double xi2 = xi*xi;
    return xi*(alpha + xi2*(beta + xi2*(gamma + xi2*delta)));
}

#include <limits>

//
// Predictable normal distribution random numbers
//
// Instead of the std::normal_distribution, which uses a rejection
// scheme, so one cannot know how many times the underlying generator
// has been called, this class requests exactly two random numbers to
// be generated for every two normal random numbers returned.
//
// This makes skip ahead parallelization reproducable
//
class predictable_normal_distribution
{
  public:
    using result_type = double;
    struct param_type {
        typedef predictable_normal_distribution distribution_type;
        param_type(): mean_{0.0}, stddev_{0.0}, precise_{false} {}
        explicit param_type(double mean, double stddev, bool precise)
        : mean_(mean), stddev_(stddev), precise_(precise)
        {}
        double mean() const {
            return mean_;
        }
        double stddev() const {
            return stddev_;
        }
        double precise() const {
            return precise_;
        }
        bool operator==(const param_type& other) const {
            return mean_ == other.mean_ && stddev_ == other.stddev_ && precise_ == other.precise_;
        }
      private:
        double mean_, stddev_;
        bool precise_;
    };     
  private:
    const param_type p_;
    bool have_x2_;
    result_type x2_;
  public:    
    explicit predictable_normal_distribution(double mean=0.0, double stddev=1.0, bool precise=false)
    : p_{mean, stddev,precise}, have_x2_{false}, x2_{}
    {}
    void reset() {
        have_x2_ = false;
        x2_   = double{};
    }
    double mean() const {
        return p_.mean();
    }
    double stddev() const {
        return p_.stddev();
    }
    double precise() const {
        return p_.precise();
    }
    param_type param() const {
        return p_;
    }
    constexpr result_type min() const {
        return std::numeric_limits<double>::lowest();
    }
    constexpr result_type max() const {
        return std::numeric_limits<double>::max();
    }
    template<class URNG>
    std::pair<result_type,result_type> draw_two(URNG& g, const param_type& params)
    {
        // generate a pair of random variables
        drand48dist dist;
        const result_type y1 = dist(g); 
        const result_type y2 = dist(g);
        const result_type fac = sqrt(-2*log(y1));
        constexpr const double pi = 3.1415926535897932;
        result_type x2 = params.precise()?cos(2*pi*y2):fastcostwopi(y2);
        result_type x1 = sqrt(1-x2*x2);
        if (y2 > 0.5)
            x1 = -x1;
        x1 *= fac;
        x2 *= fac;
        return {params.mean() + params.stddev()*x1,
                params.mean() + params.stddev()*x2};
    }
    template<class URNG>
    result_type operator()(URNG& g, const param_type& params)
    {
        static param_type  previous_params;
        static bool        local_have_x2{false};
        static result_type local_x2{0};
        if (params==previous_params && local_have_x2) { // already one available ?
            local_have_x2 = false;
            return local_x2;
        } else {
            const std::pair<result_type,result_type> x = draw_two(g, params);
            local_have_x2 = true;
            local_x2 = x.second;
            previous_params = params;
            return x.first;
        }
    }
    template<class URNG>
    result_type operator()(URNG& g)
    {
        if (have_x2_) {  // already one available ?
            have_x2_ = false;
            return x2_;
        } else {
            const std::pair<result_type,result_type> x = draw_two(g, p_);
            x2_ = x.second;
            have_x2_ = true;
            return x.first;
        }
    }
    bool operator==(const predictable_normal_distribution& other) const {
        return p_ == other.p_ && have_x2_ == other.have_x2_ && x2_ == other.x2_;
    }
    bool operator!=(const predictable_normal_distribution& other) const {
        return ! this->operator==(other);
    }
};

#endif
