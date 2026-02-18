///
/// @file debug.h
///
/// @brief TDB
///
/// @author Ramses van Zon
/// @date 2026
///

#ifndef DEBUG_H
#define DEBUG_H

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <unistd.h>

#ifndef NDEBUG
 #ifndef ATTACHCHECKUBUNTU
  #define ATTACHCHECKUBUNTU {char blockingfile[]="/proc/sys/kernel/yama/ptrace_scope"; FILE* f=fopen(blockingfile,"r"); if (f) { char c=getc(f); if (c=='1') fprintf(stderr,"sudo sh -c 'echo 0 > %s'\n",blockingfile); fclose(f);} }
 #endif
 #ifndef EMACSDEBUG
  #define EMACSDEBUG {ATTACHCHECKUBUNTU; int i=0;fprintf(stderr,"\nemacs -q -fn 6x12 -T rank%d --eval '(gdb \"gdb --annotate=1 -q -i=mi -ex finish -ex finish -ex \\\"set var i=1\\\" -ex next ./pardy -p %d\")'&\n",global_rank,getpid());while(0==i)sleep(5);}
 #endif
 #ifndef XTERMDEBUG
  #define XTERMDEBUG {ATTACHCHECKUBUNTU; int i=0;fprintf(stderr,"\nxterm -fn 6x12 -geometry 100x25 -T rank%d -e 'gdb -tui -q -ex finish -ex finish -ex \"set var i=1\" -ex next ./pardy -p %d'&\n",global_rank,getpid());while(0==i)sleep(5);}
 #endif
#else
 #ifndef EMACSDEBUG
  #define EMACSDEBUG do{;}while(0)
 #endif
 #ifndef XTERMDEBUG
  #define XTERMDEBUG do{;}while(0)
 #endif
#endif

#define ENTERDEBUGGER  \
    if (argc>1 && strncmp(argv[1],"-e", 3)==0) {                    \
        argv++; argc--; EMACSDEBUG;                       \
    } else if (argc>1 && strncmp(argv[1],"-x",3)==0) {              \
        argv++; argc--; XTERMDEBUG;                       \
    }

#ifndef NDEBUG

// specialized asserts 

#define assert_eq(a,b) if(!((long long)(a)==(long long)(b))){std::cerr << #a << ':' << (long long)(a) << ' ' << #b << ':' << (long long)(b) << '@' << global_rank << std::endl;assert((long long)(a)==(long long)(b));}
#define assert_neq(a,b) if(!((long long)(a)!=(long long)(b))){std::cerr << #a << ':' << (long long)(a) << ' ' << #b << ':' << (long long)(b) << '@' << global_rank << std::endl;assert((long long)(a)!=(long long)(b));}
#define assert_lt(a,b) if(!((long long)(a)<(long long)(b))){std::cerr << #a << ':' << (long long)(a) << ' ' << #b << ':' << (long long)(b) << '@' << global_rank << std::endl;assert((long long)(a)<(long long)(b));}
#define assert_le(a,b) if(!((long long)(a)<=(long long)(b))){std::cerr << #a << ':' << (long long)(a) << ' ' << #b << ':' << (long long)(b) << '@' << global_rank << std::endl;assert((long long)(a)<=(long long)(b));}
#define assert_gt(a,b) if(!((long long)(a)>(long long)(b))){std::cerr << #a << ':' << (long long)(a) << ' ' << #b << ':' << (long long)(b) << '@' << global_rank << std::endl;assert((long long)(a)>(long long)(b));}
#define assert_ge(a,b) if(!((long long)(a)>=(long long)(b))){std::cerr << #a << ':' << (long long)(a) << ' ' << #b << ':' << (long long)(b) << '@' << global_rank << std::endl;assert((long long)(a)>=(long long)(b));}

#else

#define assert_eq(a,b)
#define assert_neq(a,b)
#define assert_lt(a,b)
#define assert_le(a,b)
#define assert_gt(a,b)
#define assert_ge(a,b)

#endif

#endif
