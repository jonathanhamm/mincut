#ifndef _BIS_H_
#define _BIS_H_

#include <stdint.h>

#define _MUTATIONPROB 5
#define _POOLSIZE 50
#define _GET_CHRBYSIZE(pool) ((pool->chromsize >> 3) + \
                                  (pool->chromsize % 8 != 0))
#define _GET_CHBITLEN(pool) ((pool->chromsize << 6) - (64 - pool->remain))
#define _CQWORDSIZE(csize) (((csize) >> 6) + ((csize) % 64 != 0)) 

typedef struct pool_s pool_s;

/*Chromosomes are Little Endian, and Packed into a 64-bit alligned buffer*/
typedef struct roulette_s roulette_s;

struct roulette_s
{
  float prob;
  uint64_t *ptr;
};

struct pool_s
{
  uint64_t gen;       /*generation number*/
  uint16_t chromsize; /*size in quad words*/
  uint64_t cmask;     /*bit mask for each chromosome*/
  uint8_t remain;     /*Remainder bits that carry over in new quad word*/
  pool_s *parent;
  pool_s *child;
  wgraph_s *graph;
  uint64_t *crbackup;
  uint64_t *crmask;
  void (*cross) (pool_s *, uint64_t *, uint64_t *);
  void (*mutate) (pool_s *, uint64_t *);
  roulette_s rbuf[_POOLSIZE];
  uint64_t popul[0];
};


/*ge routines*/
extern void printpool (pool_s *p);
extern int run_ge (wgraph_s *g);

extern int countdigits(pool_s *p, uint64_t *cptr);
extern void printweights (pool_s *p);

/* Crossover Functions */
void singlepoint_cr (pool_s *p, uint64_t *p1, uint64_t *p2);
void mask_cr (pool_s *p, uint64_t *p1, uint64_t *p2);


/* Mutation Function */
void mutate1 (pool_s *p, uint64_t *victim);

#endif
