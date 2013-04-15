#ifndef _BIS_H_
#define _BIS_H_

#include <stdint.h>
#include <time.h>

#define _INITMUTATIONPROB 30
#define _POOLSIZE 50
#define _GET_CHRBYSIZE(pool) ((pool->chromsize / 8) + \
                                  (pool->chromsize % 8 != 0))
#define _GET_CHBITLEN(pool) (pool->bitlen)
#define _CQWORDSIZE(csize) (((csize) / 64) + ((csize) % 64 != 0)) 

typedef struct pool_s pool_s;

/*Chromosomes are Little Endian, and Packed into a 64-bit alligned buffer*/
typedef struct roulette_s roulette_s;

struct roulette_s
{
  float prob;
  float fitness;
  uint64_t *ptr;
};

struct pool_s
{
  uint64_t gen;       /*generation number*/
  uint16_t chromsize; /*size in quad words*/
  uint16_t bitlen;    /*size in bits*/
  uint64_t cmask;     /*bit mask for each chromosome*/
  uint8_t remain;     /*Remainder bits that carry over in new quad word*/
  pool_s *parent;
  pool_s *child;
  uint8_t mutateprob;
  wgraph_s *graph;
  uint64_t *crbackup;
  uint64_t *crmask;
  clock_t start;
  void (*cross) (pool_s *, roulette_s *, roulette_s *);
  void (*mutate) (pool_s *, uint64_t *);
  uint32_t accum;
  roulette_s rbuf[_POOLSIZE];
  uint64_t popul[0];
};


/*ge routines*/
extern void printpool (pool_s *p);
extern int run_ge (wgraph_s *g);

extern int countdigits(pool_s *p, uint64_t *cptr);
extern void printweights (pool_s *p);

/* Crossover Functions */
extern void singlepoint_cr (pool_s *p, uint64_t *p1, uint64_t *p2);
extern void mask_cr (pool_s *p, roulette_s *rp1, roulette_s *rp2);


/* Mutation Functions */
extern void mutate1 (pool_s *p, uint64_t *victim);
extern void mutate2 (pool_s *p, uint64_t *victim);


#endif
