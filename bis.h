#ifndef _BIS_H_
#define _BIS_H_

#include <stdint.h>

#define _POOLSIZE 50
#define _GET_CHRBYSIZE(pool) ((pool->chromsize >> 3) + \
                                  (pool->chromsize % 8 != 0))

typedef struct pool_s pool_s;

struct pool_s {
  uint32_t gen;       /*generation number*/
  uint16_t chromsize; /*size in bits*/
  uint64_t cmask;     /*bit mask for each chromosome*/
  pool_s *parent;
  pool_s *child;
  uint64_t *popul;
};

extern pool_s *pool_s_ (uint16_t csize);

extern unsigned char *read_gfile(const char *fname);

/*Graph Parsing Routines*/

/*  Graph Grammar:
 *  <graph> => <nodelist> <edgelist> EOF
 *  <nodelist> => N={n <nodeparam>}
 *  <nodeparam> => ,n <nodeparam> | E
 *  <edgelist> => E={<e> <edgeparam> }
 *  <edgeparam> => ,<e> <edgeparam> | E 
 *  <e> => {n,n,real}
 *
 *  Regex for node/edge definitions:
 *  <N> : (a...Z)+(a...z+0...9)*
 *  <e> : (0...9)+<optional_fract>
 *  <optional_frack> : (0...9)?
 */

#endif
