#ifndef DEBUG_H
#define DEBUG_H

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

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

/* Print functions used in specialized asserts; these leak memory
   like crazy, but should only be used just before the program crashes
   anyway. */

#define maxoutchars 60
static char* floatprint(double x)
{
    char* s=malloc(maxoutchars);
    snprintf(s, maxoutchars, "%e", x);
    return s;
}
static char* doubleprint(double x)
{
    char* s=malloc(maxoutchars);
    snprintf(s, maxoutchars, "%e", x);
    return s;
}
static char* intprint(int x){
    char* s=malloc(maxoutchars);
    snprintf(s, maxoutchars, "%d", x);
    return s;
}
static char* longprint(long x){
    char* s=malloc(maxoutchars);
    snprintf(s, maxoutchars, "%ld", x);
    return s;
}
static char* longunsignedprint(long unsigned x){
    char* s=malloc(maxoutchars);
    snprintf(s, maxoutchars, "%lu", x);
    return s;
}
/* static char* mpirequestprint(int n, MPI_Request x[]) { */
/*     char* s=malloc(n*maxoutchars); */
/*     char* r=s; */
/*     s += snprintf(s, 2, "{"); */
/*     for (int i=0;i<n;i++) */
/*         s += snprintf(s, maxoutchars, "[%d]=%p, ", i,x[i]); */
/*     s += snprintf(s-2, 2, "}"); */
/*     return r; */
/* } */

#define printg(x) _Generic((x), int: intprint(x), long: longprint(x), long unsigned: longunsignedprint(x), float: floatprint(x), double: doubleprint(x), char: (&(x)), char*: x)

/* specialize asserts */

#define assert_eq(a,b) if(!(a==b)){fprintf(stderr,"%s:%s %s:%s @%d\n",#a,printg(a),#b,printg(b),global_rank);assert(a==b);}
#define assert_neq(a,b) if(!(a!=b)){fprintf(stderr,"%s:%s %s:%s @%d\n",#a,printg(a),#b,printg(b),global_rank);assert(a!=b);}
#define assert_lt(a,b) if(!(a<b)){fprintf(stderr,"%s:%s %s:%s @%d\n",#a,printg(a),#b,printg(b),global_rank);assert(a<b);}
#define assert_le(a,b) if(!(a<=b)){fprintf(stderr,"%s:%s %s:%s @%d\n",#a,printg(a),#b,printg(b),global_rank);assert(a<=b);}
#define assert_gt(a,b) if(!(a>b)){fprintf(stderr,"%s:%s %s:%s @%d\n",#a,printg(a),#b,printg(b),global_rank);assert(a==b);}
#define assert_ge(a,b) if(!(a>=b)){fprintf(stderr,"%s:%s %s:%s @%d\n",#a,printg(a),#b,printg(b),global_rank);assert(a>=b);}
#else
#define assert_eq(a,b)
#define assert_neq(a,b)
#define assert_lt(a,b)
#define assert_le(a,b)
#define assert_gt(a,b)
#define assert_ge(a,b)
#endif

#endif
