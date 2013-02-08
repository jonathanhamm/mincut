#ifndef _BIS_H_
#define _BIS_H_

#include <stdint.h>

#define _POOLSIZE 50
#define _GET_CHRBYSIZE(pool) ((pool->chromsize >> 3) + \
                                  (pool->chromsize % 8 != 0))

typedef struct pool_s pool_s;

struct pool_s 
{
  uint32_t gen;       /*generation number*/
  uint16_t chromsize; /*size in bits*/
  uint64_t cmask      /*bit mask for each chromosome*/
  pool_s *parent;
  pool_s *child;
  uint64_t *popul;
};


#endif
