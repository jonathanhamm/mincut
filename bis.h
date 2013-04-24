#ifndef _BIS_H_
#define _BIS_H_

#include <stdint.h>
#include <time.h>

#define CR_N 2
#define TOURN_K 85
#define CBUF_SIZE 32
#define INITMUTATIONPROB 3
#define POOLSIZE 50
#define CRBACKUP1 POOLSIZE
#define CRBACKUP2 (POOLSIZE+1)

#define SIMA_i0     100
#define SIMA_t0     500
#define SIMA_alpha  0.80
#define SIMA_beta   1.05
#define SIMA_curr   0
#define SIMA_tmp    1
#define SIMA_RAND() (float)rand()/(float)RAND_MAX

#define SIMULATED_ANNEALING 1
#define HILL_CLIMBNG  0

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
  double prob;
  double fitness;
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
  uint8_t mutateprob;
  wgraph_s *graph;
  uint64_t start;
  void (*select) (selected_s *);
  void (*cross) (uint64_t *, uint64_t *, uint64_t *, uint64_t *);
  void (*mutate) (uint64_t *);
  int  (*e_pow) (void); /* The power calculation used in simluated annealing which is absent from hill climbing*/
  double fitsum;
  uint16_t ranksum;
  uint32_t accum;
  uint8_t k;
  uint64_t *bestfeasible;
  roulette_s rbuf[POOLSIZE];
  /* Simulated Annealing and Hill Climbing */
  float iterations, T, alpha, beta;
#define is_not_ge e_pow
#define perturb mutate
#define solusize chromsize
#define solution popul
  uint64_t popul[0];
};

extern pool_s *pool_;

/*Genetic Algorithm routines*/
extern void printpool (void);
extern int run_ge (wgraph_s *g);

extern int countdigits(uint64_t *cptr);
extern void printweights (void);

/* Selection Functions */
extern void roulette_sf (selected_s *parents);
extern void rank_sf (selected_s *parents);
extern void tournament_sf (selected_s *parents);


/* Crossover Functions */
extern void npoint_cr (uint64_t *p1, uint64_t *p2, uint64_t *dst1, uint64_t *dst2);
extern void uniform_cr (uint64_t *p1, uint64_t *p2, uint64_t *dst1, uint64_t *dst2);



/* Mutation Functions */
extern void mutate1 (uint64_t *victim);
extern void mutate2 (uint64_t *victim);

extern void printsolution (int index);
extern void printgestatus (void);
extern void printsastatus (void);

/*simulated annealing functions*/
extern int run_simanneal (wgraph_s *g, int sa_hc);

#endif
