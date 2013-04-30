/*
 bis.h
 Author: Jonathan Hamm
 
 Description:
 
 Library implemented by bis.c
 This is a library with all the functions that run the
 genetic algorithm, simulated annealing, and foolish 
 hill climbing. The code for these algorithms is 
 implemented in bis.c
 */

#ifndef _BIS_H_
#define _BIS_H_

#include <stdint.h>
#include <time.h>

/*  
    Macros soley used for inserting testing code for the algorithms,
    and termination criteria.
 */
#define MAX_GENGUESS(val) (20000 * pool_->bitlen)
#define TESTMODE
//#undef TESTMODE
#ifdef TESTMODE
    #define OPTIMAL 67
    #define MGEN
    #undef MGEN
    #define MAX_GENERATIONS MAX_GENGUESS(pool_->bitlen)
    #define EXITOPTGE()  {\
                        printf ("\nFound Optimal\n"); \
                        printgestatus (); \
                        printsolution (-1, pool_->bestfeasible); \
                        kill(getppid(), SIGQUIT); \
                        exit(EXIT_SUCCESS);\
                        }
    #define EXITGENGE()   {\
                        printf ("\nMaximum # Generations Reached\n"); \
                        printgestatus (); \
                        printsolution (-1, pool_->bestfeasible); \
                        kill(getppid(), SIGQUIT); \
                        exit(EXIT_SUCCESS);\
                        }

    #define EXITOPTSA() {\
                        printf ("\nFound Optimal\n"); \
                        printsastatus (); \
                        printsolution (-1, pool_->bestfeasible); \
                        kill(getppid(), SIGQUIT); \
                        exit(EXIT_SUCCESS);\
                        }
    #define EXITGENSA() {\
                        printf ("\nMaximum # Generations Reached\n"); \
                        printsastatus (); \
                        printsolution (-1, pool_->bestfeasible); \
                        kill(getppid(), SIGQUIT); \
                        exit(EXIT_SUCCESS);\
                        }
#endif

/* Genetic algorithm constatns and macros */
#define CR_N 5
#define MDIV_CONST 5
#define NM_BITS 2
#define TOURN_K 85
#define CBUF_SIZE 32
#define INITMUTATIONPROB 5
#define POOLSIZE 35
#define CRBACKUP1 POOLSIZE
#define CRBACKUP2 (POOLSIZE+1)
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

/* Simulated Annealing Constants */
#define SIMA_i0     500
#define SIMA_t0     4000
#define SIMA_alpha  0.95     /* 0 < alpha < 1 */
#define SIMA_beta   1.03     /* beta > 1, but preferable: 1.01 <= beta <= 1.05 */
#define SIMA_curr   0
#define SIMA_tmp    1
/* Generates a random float between 0 and 1 */
#define SIMA_RAND() (float)rand()/(float)RAND_MAX

#define SIMULATED_ANNEALING 1
#define HILL_CLIMBNG  0

typedef struct pool_s pool_s;
typedef struct roulette_s roulette_s;
typedef struct selected_s selected_s;
typedef struct ppair_s ppair_s;

/* Structure used for roullette and rank */
struct roulette_s
{
    double prob;
    double fitness;
    uint32_t cummulative;
    uint64_t *ptr;
};

/* Parent Pair (for selection)*/
struct ppair_s
{
    roulette_s *p1;
    roulette_s *p2;
};

/* Array of parents selected. */
struct selected_s
{
    ppair_s couples[NSELECT];
};

/* Main data structure for genetic algorithm, simulated annealing, and foolish hill climbing */
struct pool_s
{
    uint64_t gen;               /*generation number*/
    uint16_t nqwords;           /*size of chromosome in quad (64 bit) words*/
    uint16_t bitlen;            /*size of chromosome in bits*/
    uint64_t cmask;             /*bit mask for each chromosome*/
    uint8_t remain;             /*Remainder bits that carry over in new quad word*/
    uint8_t mutateprob;         /*mutation probability*/
    wgraph_s *graph;            /*Graph data structure that this references*/
    uint64_t start;             /*start time*/
    double fitsum;              /*sum of all fitnesses*/
    uint16_t ranksum;           /*sum of all ranks*/
    uint32_t accum;             /*sum of selection probabilities*/
    uint8_t k;                  /*k value used in tournament selection.*/
    uint64_t *bestfeasible;     /*Pointer to best feasible Chromosome*/
    roulette_s rbuf[POOLSIZE];  /*Sorted array (by probability) of Chromosome pointers*/
    
    /* Dynamically Bound Functions: */
    
    /* Simulated Annealing and Hill Climbing Parameters*/
    double iterations, T, alpha, beta;
    /* Selection Function */
    void (*select) (selected_s *);
    /* Crossover Function */
    void (*cross) (uint64_t *, uint64_t *, uint64_t *, uint64_t *);
    /* Mutation Function */
    void (*mutate) (uint64_t *);
     /* 
      This is a 'boolean' function that represents the boolean expression OR'd with
      with the comparrison of the newly generated solution in simulated annealing: 
        if( h(NewS) < h(S) or random < e ^ [ (h(S) - h(NewS)) / T ] )
        -This function is a boolean for the right side of the OR statement.
      For simulated annealing, this function contains the power expression. For foolish
      hill climbing, this is overloaded with a NOP function that always returns 0 (false). 
      */
    int  (*e_pow) (void);
#define is_not_ge e_pow
#define perturb mutate
#define solusize nqwords
#define solution popul
#define nperturbations gen
    /* Population. This is contiguously appended at the end of the structure at runtime. */
    uint64_t popul[0];
};

/* Global Pool */
extern pool_s *pool_;

/*Genetic Algorithm routines*/
extern void printpool (void);
extern int run_ge (wgraph_s *g);

/* Selection Functions */
extern void roulette_sf (selected_s *parents);
extern void rank_sf (selected_s *parents);
extern void tournament_sf (selected_s *parents);

/* Crossover Functions */
extern void npoint_cr (uint64_t *p1, uint64_t *p2, uint64_t *dst1, uint64_t *dst2);
extern void uniform_cr (uint64_t *p1, uint64_t *p2, uint64_t *dst1, uint64_t *dst2);

/* Mutation/Perturbation Functions */
extern void mutate1 (uint64_t *victim);
extern void mutate2 (uint64_t *victim);
extern void pairwise_ex (uint64_t *victim);
extern void perturbinvert (uint64_t *victim);


/* Printing Functions */
extern void printsolution (int index, uint64_t *ptr);
extern void printgestatus (void);
extern void printsastatus (void);

/*simulated annealing functions*/
extern int run_simanneal (wgraph_s *g, int sa_hc);

#endif