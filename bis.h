#ifndef _BIS_H_
#define _BIS_H_

#include <stdint.h>

#define _POOLSIZE 50
#define _GET_CHRBYSIZE(pool) ((pool->chromsize >> 3) + \
                                  (pool->chromsize % 8 != 0))

#define _CQWORDSIZE(csize) (((csize) >> 6) + ((csize) % 64 != 0)) 

typedef struct pool_s pool_s;

/*Chromosomes are Little Endian, and Packed into a 64-bit alligned buffer*/
typedef struct chromosome_s chromosome_s;

struct pool_s
{
  uint32_t gen;       /*generation number*/
  uint16_t chromsize; /*size in bits*/
  uint64_t cmask;     /*bit mask for each chromosome*/
  pool_s *parent;
  pool_s *child;
  uint64_t popul[0];
};

pool_s *pool_init (wgraph_s *g);

/*ge routines*/
extern pool_s *pool_init (wgraph_s *g);
extern void printpool (pool_s *p);

#endif
