#ifndef _BIS_H_
#define _BIS_H_

#include <stdint.h>
#include <time.h>

#define CBUF_SIZE 32
#define INITMUTATIONPROB 3
#define POOLSIZE 50
#define PSIZE_ROUL_DIV 5

#if POOLSIZE / PSIZE_ROUL_DIV == 0
  #define NSELECT 1
#else
  #define NSELECT (POOLSIZE / PSIZE_ROUL_DIV)
#endif

#define GET_CHRBYSIZE(pool) ((pool->chromsize / 8) + \
                                  (pool->chromsize % 8 != 0))
#define GET_CHBITLEN(pool) (pool->bitlen)
#define CQWORDSIZE(csize) (((csize) / 64) + ((csize) % 64 != 0)) 

typedef struct pool_s pool_s;

/*Chromosomes are Little Endian, and Packed into a 64-bit alligned buffer*/
typedef struct roulette_s roulette_s;
typedef struct selected_s selected_s;
typedef struct ppair_s ppair_s;


struct roulette_s
{
  float prob;
  float fitness;
  uint32_t cummulative;
  uint64_t *ptr;
};

struct ppair_s
{
  roulette_s *p1;
  roulette_s *p2;
};

struct selected_s
{
  ppair_s couples[NSELECT];
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
  uint64_t start;
  void (*select) (pool_s *, selected_s *);
  void (*cross) (pool_s *, roulette_s *, roulette_s *);
  void (*mutate) (pool_s *, uint64_t *);
  float fitsum;
  uint32_t accum;
  uint64_t *bestfeasible;
  roulette_s rbuf[POOLSIZE];
  uint64_t popul[0];
};


/*ge routines*/
extern void printpool (pool_s *p);
extern int run_ge (wgraph_s *g);

extern int countdigits(pool_s *p, uint64_t *cptr);
extern void printweights (pool_s *p);

/* Selection Functions */
extern void roulette_sf (pool_s *p, selected_s *parents);

/* Crossover Functions */
extern void singlepoint_cr (pool_s *p, uint64_t *p1, uint64_t *p2);
extern void uniform_cr (pool_s *p, roulette_s *rp1, roulette_s *rp2);


/* Mutation Functions */
extern void mutate1 (pool_s *p, uint64_t *victim);
extern void mutate2 (pool_s *p, uint64_t *victim);

extern void buildgraph (pool_s *p);

#endif
