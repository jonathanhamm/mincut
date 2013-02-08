#ifndef _BIS_H_
#define _BIS_H_

#include <stdint.h>

#define GET_CHROMBUFFER(pool) ((uint64_t *)(pool->child + sizeof(pool_s *)))
#define GET_CHROMBYTESIZE(pool) ((pool->chromsize >> 3) + \
                                  (pool->chromsize % 8 != 0))

typedef struct pool_s pool_s;

struct pool_s 
{
  uint32_t gen;       /*generation number*/
  uint16_t chromsize; /*size in bits*/
  uint16_t popsize;   /*number of chromosomes*/
  uint64_t cmask      /*bit mask for each chromosome*/
  pool_s *parent;
  pool_s *child;
  /*Allocated at Runtime:*/
  /*uint64_t chromosomes[...]*/ 
};


#endif
