/*
 bis.c
 Author: Jonathan Hamm
 
 Description: 
 
 This file contains the main genetic algorithm code, simulated 
 annealing code, and foolish hill climbing code. This includes 
 selection functions, crossover functions, mutation operators, 
 fitness functions, perturbation functions, and the main loops 
 for the genetic algorithm, simulated annealing, and foolish
 hill climbing. 
 */

#include "parse.h"
#include "bis.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <unistd.h>
#include <sys/sem.h>
#include <signal.h>

/*
 Counts the set bits in a 32-bit long word:
 32-bit count code obtained from:
 http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
 */
#define countbitsLW(lw) ((lw & 0xfff) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f \
+ (((lw & 0xfff000) >> 12) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f \
+ ((lw >> 24) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f

typedef struct solset_s solset_s;

struct solset_s
{
    vertex_s *v;
    solset_s *next;
};

/* Global Variables */
static pid_t pid_;
static int pipe_[2];
pool_s *pool_;

/* Constructors for Genetic Algorithm Pool and Simulated Annealing & Hill Climbing "Pool" */
static void pool_s_ (wgraph_s *g);
static void pool_s_simanneal (wgraph_s *g, int sa_hc);

/* Bit Operations */
static uint8_t getbit (uint64_t *chrom, uint16_t pos);
static void setbit (uint64_t *chrom, uint16_t pos, uint8_t val);

/* Fitness Evaluation Functions */
static uint16_t countdigits(uint64_t *cptr);
static int iscut(uint64_t *chrom, vertex_s *v);
static double sumweights (uint64_t *chrom);
static double getfitness (uint64_t *chrom);
static int prcmp (roulette_s *a, roulette_s *b);
static void computeprob (void);
static int isfeasible (uint64_t *chrom);

/* Iterative Binary Search */
static int bsearch_r (roulette_s *roul, uint32_t key);

/* Simulated Annealing and Foolish Hill Climbing Functions */
static void sima_rand (roulette_s *dst);
static int e_pow_sa (void);
static int nop_hc (void);

/* Signal Handlers */
static void cSIGUSR1 (int signal);
static void pSIGINT (int signal);
static void sigNOP (int signal);

/* Printing Functions */
static void printqword (uint64_t lword, uint8_t mask);
static void printchrom (uint64_t *chrom);
static void insert_solset (solset_s **sset, vertex_s *v);
static void print_solset (solset_s *solset);

/*
 "Constructor" for pool_s structure. This simply initializes a pool
 of chromosomes for the genetic algorithm.
 
 @param g   Pointer to graph data structure used to initialize pool.
 */
void pool_s_ (wgraph_s *g)
{
    uint16_t  i, j;
    uint64_t  *ptr;
    double sum;
    
    pool_ = calloc(1, sizeof(*pool_) + POOLSIZE * CQWORDSIZE(g->nvert+1) * 8);
    if (!pool_) {
        perror ("Heap Allocation Error");
        exit (EXIT_FAILURE);
    }
    pool_->chromsize = (g->nvert / 64) + (g->nvert % 64 != 0);
    pool_->remain = g->nvert % 64;
    /* 
     For chromosomes that don't align with 64 bits, the extra bits
     in the partially filled (last quad word) 64-bit word need to be
     masked out. This code creates the mask. 
     */
    if (!(g->nvert % 64) && g->nvert)
        pool_->cmask = 0xffffffffffffffffllu;
    else {
        for (i = 0; i < pool_->remain; i++)
            pool_->cmask |= (1llu << i);
    }
    pool_->select = tournament_sf;
    pool_->cross = uniform_cr;
    pool_->mutate = mutate1;
    pool_->k = TOURN_K;
    pool_->gen = 0;
    ptr = pool_->popul;
    /* Initialize Chromosomes Randomly */
    for (i = 0; i < POOLSIZE; i++) {
        for (j = 0; j < pool_->chromsize; j++, ptr++) {
            ((uint16_t *)ptr)[0] = (uint16_t)rand();
            ((uint16_t *)ptr)[1] = (uint16_t)rand();
            ((uint16_t *)ptr)[2] = (uint16_t)rand();
            ((uint16_t *)ptr)[3] = (uint16_t)rand();
        }
        *(ptr-1) &= pool_->cmask;
    }
    pool_->graph = g;
    pool_->mutateprob = INITMUTATIONPROB;
    pool_->bitlen = ((pool_->remain) ? ((pool_->chromsize * 64) - (64 - pool_->remain)) : pool_->chromsize*64);
    for (sum = 0, ptr = pool_->popul, i = 0; i < POOLSIZE; i++, ptr += pool_->chromsize) {
        pool_->rbuf[i].fitness = getfitness (ptr);
        pool_->rbuf[i].ptr = ptr;
        sum += pool_->rbuf[i].fitness;
    }
    for (i = 0, pool_->accum = 0; i < POOLSIZE; i++)
        pool_->rbuf[i].prob = sum / pool_->rbuf[i].fitness;
    qsort (pool_->rbuf, POOLSIZE, sizeof(roulette_s), (int (*)(const void *, const void *))prcmp);
    for (i = 0; i < POOLSIZE; i++) {
        pool_->accum += (int)pool_->rbuf[i].prob;
        pool_->rbuf[i].cummulative = (i+1);
        if (i)
            pool_->rbuf[i].cummulative += pool_->rbuf[i-1].cummulative;
    }
    pool_->ranksum = pool_->rbuf[POOLSIZE-1].cummulative;
    pool_->fitsum = sum;
}

/*
 "Constructor" for simulated annealing and hill climbing "pool" (1 chromosome). 
 
 @param g       Pointer to graph data structure used to initialize pool.
 @param sa_hc   A boolean indicating if this is simulated annealing or hill climbing. 
 */
void pool_s_simanneal (wgraph_s *g, int sa_hc)
{
    int i;
    
    pool_ = calloc (1, sizeof(*pool_) + 2 * CQWORDSIZE(g->nvert+1) * 8);
    if (!pool_) {
        perror ("Heap Allocation Error");
        exit (EXIT_FAILURE);
    }
    pool_->solusize = (g->nvert / 64) + (g->nvert % 64 != 0);
    pool_->remain = g->nvert % 64;
    pool_->bitlen = ((pool_->remain) ? ((pool_->chromsize * 64) - (64 - pool_->remain)) : pool_->chromsize*64);
    /*
     For solutions that don't align with 64 bits, the extra bits
     in the partially filled (last quad word) 64-bit word need to be
     masked out. This code creates the mask. 
     */
    if (!(g->nvert % 64) && g->nvert)
        pool_->cmask = 0xffffffffffffffffllu;
    else {
        for (i = 0; i < pool_->remain; i++)
            pool_->cmask |= (1llu << i);
    }
    pool_->perturb = mutate1;
    pool_->graph = g;
    pool_->T = SIMA_t0;
    pool_->iterations = SIMA_i0;
    pool_->alpha = SIMA_alpha;
    pool_->beta = SIMA_beta;
    if (sa_hc == SIMULATED_ANNEALING)
        pool_->e_pow = e_pow_sa;
    else
        pool_->e_pow = nop_hc;
    pool_->bestfeasible = calloc (sizeof(uint64_t), pool_->solusize);
    if (!pool_->bestfeasible) {
        perror ("Heap Allocation Error");
        exit (EXIT_FAILURE);
    }
}

/*
 Returns a bit value from a chromosome at the
 position specified by 'pos'.
 
 @param chrom   Chromosome to look at. 
 @param pos     Position in chromosome. 
 @return        Returns the bit value at position 'pos'. 
 */
uint8_t getbit (uint64_t *chrom, uint16_t pos)
{
    return (chrom[pos / 64] >> (pos % 64)) & 1llu;
}

/*
 Sets a bit from a chromosome at the position specified 
 by 'pos'.
 
 @param chrom   Chromosome to to set bit in. 
 @param pos     Position in chromosome to set bit. 
 @param val     Value to set chromosome to. 
 */
void setbit (uint64_t *chrom, uint16_t pos, uint8_t val)
{
    chrom[pos / 64] &= ~(1llu << (pos % 64));
    if (val)
        chrom[pos / 64] |= (1llu << (pos % 64));
}

/*
 Counts the number of set bits in a chromosome. 
 
 @param cptr    Chromosome to count bits for. 
 @return        Returns the number of set bits in chromosome. 
 */
uint16_t countdigits (uint64_t *cptr)
{
    int i, count, n;
    
    n = pool_->chromsize-1;
    for (i = 0, count = 0; i < n; i++) {
        count += countbitsLW((uint32_t)cptr[i]);
        count += countbitsLW((uint32_t)(cptr[i] >> 32));
    }
    count += countbitsLW((uint32_t)(cptr[i] & pool_->cmask));
    count += countbitsLW((uint32_t)((cptr[i] & pool_->cmask) >> 32));
    return count;
}

/*
 Determines if a vertex is in the set of cut vertices (in
 other words, the vertices on the 'other side' of the 
 graph bisection). If it is on the 'other side', there is
 a cut. This locates the index of the specified vertex in
 the graph data structure (graph indices correspond to the 
 bit position of the vertex in the chromosome). The index is
 obtained from a hash to alleviate search overhead. Then the 
 function checks if that index/bit position is unset in the 
 chromosome. If it is unset, there is a cut. 
 
 @param chrom   Chromosome to function is looking in. 
 @param v       Vertex function is checking. 
 @return        Returns 1 if there is a cut, and 0 if no
                cut.
 */
int iscut(uint64_t *chrom, vertex_s *v)
{
    uint16_t i, nvert;
    vertex_s **vtable;
    
    nvert = pool_->graph->nvert;
    vtable = pool_->graph->vtable;
    i = vgetindex (v);
    if ((~chrom[i / 64] & (1llu << (i % 64))))
        return 1;
    return 0;
}

/*
 Sums the weights of all cut edges in a chromosome. 
 
 Implementation Description: This works by taking 
 each quad word (64 bit word) of the chromosome, and 
 looping through all the set bits, skipping the unset
 bits. Unset bits can be skipped using a processor's 
 find-first-set-bit instruction (many architectures have
 an instruction that can find the first set bit in a 
 register, Intel's is 'bsf'). The processed set bits 
 are then masked out to zero in the loop. The call to ffsl 
 returns the position of the first set bit. The compiler 
 inlines this function so it's just the instruction. The bit 
 position can be used to access a vertex in the graph. Then
 the function looks at every edge connected to the graph, and
 checks if it's cut, summing the cut edges. 
 
 @param chrom   Chromosome function is summing the 
                weights of cut edges for. 
 @return        Returns the sum of all cut weights. 
 */
double sumweights (uint64_t *chrom)
{
    uint8_t pos;
    uint16_t i, j, csize, nedges;
    uint64_t iter,
    *ptr;
    vertex_s *v;
    edge_s **edges;
    double weight;
    
    for (weight = 0, i = 0, ptr = chrom; i < pool_->chromsize; i++, ptr++) {
        if (i == pool_->chromsize-1 && pool_->remain)
            csize = pool_->remain;
        else
            csize = 64;
        for (iter = *ptr, pos = 0; pos < csize; iter &= ~(1llu << pos)) {
            pos = ffsl(iter);
            if (!pos)
                break;
            --pos;
            v = pool_->graph->vtable[/* i*64 */ (i << 6) + pos];
            nedges = v->nedges;
            edges = v->edges;
            for (j = 0; j < nedges; j++) {
                if (iscut(chrom, (edges[j]->v1 == v) ? edges[j]->v2 : edges[j]->v1))
                    weight += edges[j]->weight;
            }
        }
    }
    return weight;
}

/*
 Gets the fitness of a chromosome. This sums the weights of all
 cut edges, and then applies a fitness penalty if the chromosome 
 isn't feasible. 
 
 @param chrom   Chromosome to get the fitness of. 
 @return        Returns the fitness of the chromosome. 
 */
double getfitness (uint64_t *chrom)
{
    int setcount, differ;
    
    setcount = countdigits (chrom);
    differ = abs(2*setcount - pool_->bitlen);
    return sumweights (chrom) + (differ << 4);
}

/*
 This is just a callback function that is passed to
 stdlib's qsort (quicksort) function. This compares
 the probabilities of chromosomes being selected. 
 
 @param a   First structure containing a chromosome and
            its probability. 
 @param b   Second structure containing a chromosome and
            its probability. 
 @return    Returns -1 if a < b, 1 if a > b, and 0 if they
            are equal.
 */
int prcmp (roulette_s *a, roulette_s *b)
{
    if (a->prob < b->prob)
        return -1;
    if (a->prob > b->prob)
        return 1;
    return 0;
}

/*
 Computes the probability of all chromosome of being selected. Sorts
 chromosomes based on probability. 
 */
void computeprob (void)
{
    uint16_t i, n;
    double sum;
    
    sum = pool_->fitsum;
    for (pool_->accum = 0, i = 0; i < POOLSIZE; i++) {
        pool_->rbuf[i].prob = sum / pool_->rbuf[i].fitness;
        pool_->accum += (int)pool_->rbuf[i].prob;
    }
    qsort (pool_->rbuf, POOLSIZE, sizeof(roulette_s), (int (*)(const void *, const void *))prcmp);
    if (isfeasible(pool_->rbuf[POOLSIZE-1].ptr))
        pool_->bestfeasible = pool_->rbuf[POOLSIZE-1].ptr;
}

/* 
 Checks of a chromosome is feasible or not. 
 
 @param chrom   Chromosome to check if feasible. 
 @return        Returns 1 if feasible, and 0 if not. 
 */
int isfeasible (uint64_t *chrom)
{
    if (pool_->bitlen % 2)
        return  (2 * countdigits (chrom) == pool_->bitlen - 1)
        ||
        (2 * countdigits (chrom) == pool_->bitlen + 1);
    else
        return  (2 * countdigits (chrom) == pool_->bitlen);
}

/*
 Roulette selection function. 
 
 @param parents A pointer to an array of parents
                that this method will fill when
                selecting parents. 
 */
void roulette_sf (selected_s *parents)
{
    int i, j;
    uint32_t i1, i2;
    float sum;
    
    for (i = 0; i < NSELECT; i++) {
        i1 = rand() % pool_->accum;
        i2 = rand() % pool_->accum;
        for (j = 0, sum = 0; j < POOLSIZE && sum <= (float)i1; sum += pool_->rbuf[j].prob, j++);
        i1 = (j == POOLSIZE) ? j-1 : j;
        for (j = 0, sum = 0; j < POOLSIZE && sum <= (float)i2; sum += pool_->rbuf[j].prob, j++);
        i2 = (j == POOLSIZE) ? j-1 : j;
        if (i2 == i1)
            i2 = (i2 - 1) % POOLSIZE;
        parents->couples[i].p1 = &pool_->rbuf[i1];
        parents->couples[i].p2 = &pool_->rbuf[i2];
    }
}

/*
 Iterative Binary search that returns the index of rank 'region'
 that 'key' fits into. This is used for rank selection, which indexes 
 a cumulative probabily distribution. The binary search alleviates some
 of the overhead of matching the random value generated by rank_sf. 
 
 @param roul    Array to search in. 
 @param key     Search key. 
 @return        Returns the index that will contain a selected parent.  
 */
int bsearch_r (roulette_s *roul, uint32_t key)
{
	int mid, low, high;
    
    low = 0;
    high = POOLSIZE-1;
    for (mid = low+(high-low)/2; high >= low; mid = low+(high-low)/2) {
        if (key < roul[mid].cummulative)
            high = mid - 1;
        else if (key > roul[mid].cummulative)
            low = mid + 1;
        else
            return mid;
    }
    if (mid == POOLSIZE)
        --mid;
    return mid;
}

/*
 Rank Selection function. 
 
 @param parents A pointer to an array of parents
                that this method will fill when
                selecting parents.
 */
void rank_sf (selected_s *parents)
{
    int i;
    uint32_t i1, i2;
    for (i = 0; i < NSELECT; i++) {
        i1 = bsearch_r (pool_->rbuf, rand() % pool_->ranksum);
        i2 = bsearch_r (pool_->rbuf, rand() % pool_->ranksum);
        if (i2 == i1)
            i2 = (i2 - 1) % POOLSIZE;
        parents->couples[i].p1 = &pool_->rbuf[i1];
        parents->couples[i].p2 = &pool_->rbuf[i2];
    }
}

/*
 Tournament selection function. 
 
 @param parents A pointer to an array of parents
                that this method will fill when
                selecting parents.
 */
void tournament_sf (selected_s *parents)
{
    int       i;
    uint8_t   k, R;
    uint32_t  i1a, i1b,
    i2a, i2b;
    
    for (i = 0, k = pool_->k; i < NSELECT; i++) {
        i1a = rand() % POOLSIZE;
        while ((i1b = rand() % POOLSIZE) == i1a);
        i2a = rand() % POOLSIZE;
        while ((i2b = rand() % POOLSIZE) == i2a);
        R  = rand() % 100;
        if (R < k) {
            if (i1a > i1b)
                parents->couples[i].p1 = &pool_->rbuf[i1a];
            else
                parents->couples[i].p1 = &pool_->rbuf[i1b];
            if (i2a > i2b)
                parents->couples[i].p2 = &pool_->rbuf[i2a];
            else
                parents->couples[i].p2 = &pool_->rbuf[i2b];
        }
        else {
            if (i1a < i1b)
                parents->couples[i].p1 = &pool_->rbuf[i1a];
            else
                parents->couples[i].p1 = &pool_->rbuf[i1b];
            if (i2a < i2b)
                parents->couples[i].p2 = &pool_->rbuf[i2a];
            else
                parents->couples[i].p2 = &pool_->rbuf[i2b];
        }
    }
}

/*
 Performs an n-point crossover, where 'n' is specified by the macro
 constant CR_N. 
 
 Note:  dst1 and dst2 will usually equal p1 and p2, correspdoningly,
        but these sometimes have to change in order to guaruntee 
        elitism. 
 
 @param p1      Parent 1 chromosome.
 @param p2      Parent 2 chromosome.
 @param dst1    Destination chromosome 1. 
 @param dst2    Destination chromosome 2. 
 */
void npoint_cr (uint64_t *p1, uint64_t *p2, uint64_t *dst1, uint64_t *dst2)
{
    int i, j, pindex;
    uint16_t point, inv;
    uint64_t mask;
    uint64_t *backup1, *backup2;
    
    backup1 = &pool_->popul[CRBACKUP1];
    backup2 = &pool_->popul[CRBACKUP2];
    for (i = 0; i < pool_->chromsize; i++)
        backup1[i] = p1[i];
    for (i = 0; i < pool_->chromsize; i++)
        backup2[i] = p2[i];
    point = rand() % (pool_->bitlen / CR_N);
    inv = pool_->bitlen - point;
    for (i = 0; i < CR_N; i++) {
        for (j = i * (pool_->bitlen / CR_N); j < inv; j++) {
            setbit(dst1, j, getbit(backup2, point+j));
            setbit(dst2, j, getbit(backup1, point+j));
        }
        for (j = i * (pool_->bitlen / CR_N); j < point; j++) {
            setbit(dst1, inv+j, getbit(backup2, j));
            setbit(dst2, inv+j, getbit(backup1, j));
        }
    }
}

/*
 Uniform Crossover: Performs a uniform crossover where 'mask' is the 
 randomly generated bit string used to select from parents. Instead of
 looping through each allele, this can be performed 64 bits/alleles at 
 a time (with a 64-bit processor). The truth table below shows how a 
 short expression/function was derived for this.  
 
 p1 p2  mask  | Child
 ----------------
 0  0   0     | 0
 0  0   1     | 0
 0  1   0     | 0
 0  1   1     | 1
 1  0   0     | 1
 1  0   1     | 0
 1  1   0     | 1
 1  1   1     | 1
  
 Child1 = (~mask & p1) | (mask  & p2)
 Child2 = (mask  & p1) | (~mask & p2)
 
 Note:  dst1 and dst2 will usually equal p1 and p2, correspdoningly,
 but these sometimes have to change in order to guaruntee
 elitism.
 
 @param p1      Parent 1 chromosome.
 @param p2      Parent 2 chromosome.
 @param dst1    Destination chromosome 1.
 @param dst2    Destination chromosome 2.
 */
void uniform_cr (uint64_t *p1, uint64_t *p2, uint64_t *dst1, uint64_t *dst2)
{
    uint16_t i;
    uint64_t  mask,  backup;
    
    for (i = 0; i < pool_->chromsize; i++) {
        ((uint16_t *)&mask)[0] = (uint16_t)rand();
        ((uint16_t *)&mask)[1] = (uint16_t)rand();
        ((uint16_t *)&mask)[2] = (uint16_t)rand();
        ((uint16_t *)&mask)[3] = (uint16_t)rand();
        backup = p1[i];
        dst1[i] = (~mask & backup) | (mask & p2[i]);
        dst2[i] = (mask & backup) | (~mask & p2[i]);
    }
}

/*
 First mutation operator. Randomly inverts n bits, where 'n'
 is a random number of bits out of the chromosome. 'n' is bounded
 by dividing the chromosome length by a constant 'MDIV_CONST'. The
 position of these bits are randomly generated. 
 
 @param victim  Target chromosome to be mutated. 
 */
void mutate1 (uint64_t *victim)
{
    int i, n;
    uint16_t index;
    
    n = rand() % (pool_->bitlen / MDIV_CONST);
    if (!(n % 2) && !isfeasible(victim))
        ++n;
    for (i = 0; i < n; i++) {
        index = rand() % pool_->bitlen;
        victim[index / 64] ^= (1llu << (index % 64));
    }
}

/*
 Second mutation operator. Inverts a set number of bits
 determined by the NM_BITS macro constant. The positions 
 of these bits are randomly generated. 
 
 @param victim  Target chromosome to be mutated.
 */
void mutate2 (uint64_t *victim)
{
    uint16_t index, i;
    
    for (i = 0; i < NM_BITS; i++) {
        index = rand() % pool_->bitlen;
        victim[index / 64] ^= (1llu << (index % 64));
    }
}

/*
 The Pairwise Exchange perturbation function (could also be
 used as a mutation operator). This simply swaps two bits 
 from two random positions. 
 
 @param victim  Target chromosome to be mutated/perturbed.
 */
void pairwise_ex (uint64_t *victim)
{
    int backup;
    uint16_t index1, index2;
    
    index1 = rand() % pool_->bitlen;
    while ((index2 = rand() % pool_->bitlen) == index1);
    backup = getbit(victim, index1);
    setbit(victim, index1, getbit(victim, index2));
    setbit(victim, index2, backup);
}

/*
 Main function for running genetic algorithm. The genetic algorithm 
 is forked off and runs in a different thread. The parent process
 simply sends terminal commands to the genetic algorithm via signals.
 
 @param g   Graph that the genetic algorithm finds bisection on. 
 @return    Returns a 0 on successful execution. 
 */
int run_ge (wgraph_s *g)
{
    int i;
    uint16_t n, index;
    selected_s parents;
    uint64_t *p1, *p2;
    uint64_t *dst1, *dst2;
    roulette_s *rp1, *rp2, *rpt;
    unsigned char cbuf[CBUF_SIZE];
    
    if (signal(SIGINT, pSIGINT) == SIG_ERR)
        throw_exception();
    if (signal(SIGUSR1, sigNOP) == SIG_ERR)
        throw_exception();
    pool_s_ (g);
    if (!pool_)
        throw_exception();
    if (pipe (pipe_) == -1)
        throw_exception();
    pid_ = fork();
    if (pid_ < 0)
        throw_exception();
    if (pid_) {
        /* Parent Process: This only sends commands to genetic algorithm process */
        printf ("Now running Genetic Algorithm with chromosome size: %d\n> ", pool_->bitlen);
        memset (cbuf, 0, sizeof(cbuf));
        cbuf[CBUF_SIZE-1] = UEOF;
        while (1) {
            index = 0;
            while ((cbuf[index] = (char)getchar()) != '\n') {
                if (index < CBUF_SIZE-1)
                    ++index;
                else
                    printf ("Command Length %d Exceeded\n", CBUF_SIZE-1);
            }
            cbuf[index] = UEOF;
            write(pipe_[1], cbuf, CBUF_SIZE);
            kill (pid_, SIGUSR1);
            pause();
        }
    }
    else {
        if (signal(SIGUSR1, cSIGUSR1) == SIG_ERR)
            throw_exception();
        if (signal(SIGINT, sigNOP) == SIG_ERR)
            throw_exception();
        n = pool_->chromsize;
        time((time_t *)&pool_->start);
        /* 
         *
         *   Main Genetic Algorithm Loop 
         *
         */
        while (1) {
            pool_->select (&parents);
            for (i = 0; i < NSELECT; i++) {
                rp1 = parents.couples[i].p1;
                rp2 = parents.couples[i].p2;
                p1  = rp1->ptr;
                p2  = rp2->ptr;
                if (rp1->ptr == pool_->bestfeasible) {
                    rp1 = &pool_->rbuf[0];
                    dst1 = rp1->ptr;
                }
                else
                    dst1 = p1;
                if (rp2->ptr == pool_->bestfeasible) {
                    rp2 = &pool_->rbuf[0];
                    dst2 = rp2->ptr;
                }
                else
                    dst2 = p2;
                if (dst1 == dst2) {
                    for (rpt = &pool_->rbuf[rand() % POOLSIZE];
                         rpt->ptr == pool_->bestfeasible ||
                         rpt->ptr == dst1 ||
                         rpt->ptr == dst2
                         ;rpt = &pool_->rbuf[rand() % POOLSIZE]);
                    rp2 = rpt;
                    dst2 = rpt->ptr;
                }
                pool_->fitsum -= (rp1->fitness + rp2->fitness);
                pool_->cross (p1, p2, dst1, dst2);
                if (rand() % 100 < pool_->mutateprob) {
                    pool_->mutate (dst1);
                    pool_->mutate (dst2);
                }
                rp1->fitness = getfitness (dst1);
                rp2->fitness = getfitness (dst2);
                pool_->fitsum += (rp1->fitness + rp2->fitness);
                
            }
            computeprob ();
            pool_->gen++;
#ifdef TESTMODE
            if (pool_->rbuf[POOLSIZE-1].fitness == OPTIMAL) {
                printf ("\nFound Optimal\n");
                printgestatus ();
                kill(getppid(), SIGQUIT);
                exit(EXIT_SUCCESS);
            }
#endif
        }
    }
    return 0;
    
exception_:
    perror("Program Failed Initialization");
    exit (EXIT_FAILURE);
}

/*
 Main function for running simulated annealing or foolish hill climbing. The 
 simulated annealing/foolish hill climbing loop is forked off and runs in a 
 separate process. The parent process simply sends terminal commands to the 
 algorithm via signals.
 
 @param g       Graph that the simulated annealing or foolish hill climbing
                algorithm finds bisection on.
 @param sa_hc   A "boolean" that determines whether this should run as simulated
                annealing or foolish hill climbing. 
 @return        Returns a 0 on successful execution.
 */
int run_simanneal (wgraph_s *g, int sa_hc)
{
    float i;
    int j;
    uint16_t index, ssize;
    unsigned char cbuf[CBUF_SIZE];
    roulette_s *s, *new_s;
    
    pool_s_simanneal (g, sa_hc);
    s = &pool_->rbuf[SIMA_curr];
    new_s = &pool_->rbuf[SIMA_tmp];
    s->ptr = pool_->solution;
    new_s->ptr = &pool_->solution[pool_->solusize];
    sima_rand(s);
    for (j = 0; j < pool_->solusize; j++)
        pool_->bestfeasible[j] = s->ptr[j];
    if (signal(SIGINT, pSIGINT) == SIG_ERR)
        throw_exception();
    if (signal(SIGUSR1, sigNOP) == SIG_ERR)
        throw_exception();
    if (pipe (pipe_) == -1)
        throw_exception();
    pid_ = fork();
    if (pid_ < 0)
        throw_exception();
    if (pid_) {
        /* Parent Process: This only sends commands to simulated annealing/foolish hill climbing process */
        index = 0;
        if (sa_hc == SIMULATED_ANNEALING)
            printf ("Now running Simulated Annealing with solution size: %d\n> ", pool_->bitlen);
        else
            printf ("Now running \"Foolish\" Hill Climbing with solution size: %d\n> ", pool_->bitlen);
        memset (cbuf, 0, sizeof(cbuf));
        cbuf[CBUF_SIZE-1] = UEOF;
        while (1) {
            index = 0;
            while ((cbuf[index] = (char)getchar()) != '\n') {
                if (index < CBUF_SIZE-1)
                    ++index;
                else
                    printf ("Command Length %d Exceeded\n", CBUF_SIZE-1);
            }
            cbuf[index] = UEOF;
            write(pipe_[1], cbuf, CBUF_SIZE);
            kill (pid_, SIGUSR1);
            pause();
        }
    }
    else {
        if (signal(SIGUSR1, cSIGUSR1) == SIG_ERR)
            throw_exception();
        if (signal(SIGINT, sigNOP) == SIG_ERR)
            throw_exception();
        ssize = pool_->solusize;
        time((time_t *)&pool_->start);
        /*
         *
         *   Main Simulated Annealing/Foolish Hill Climbing Loop. 
         *
         */
        while (1) {
            for (i = 0; i < pool_->iterations; i++) {
                pool_->perturb (new_s->ptr);
                new_s->fitness = getfitness(new_s->ptr);
                if (new_s->fitness < s->fitness || pool_->e_pow()) {
                    for (j = 0; j < ssize; j++)
                        s->ptr[j] = new_s->ptr[j];
                    if (isfeasible(s->ptr)) {
                        for (j = 0; j < ssize; j++)
                            pool_->bestfeasible[j] = s->ptr[j];
                    }
                    s->fitness = new_s->fitness;
                }
                pool_->nperturbations++;
            }
            pool_->T = pool_->alpha * pool_->T;
            pool_->iterations = pool_->beta * pool_->iterations;
        }
    }
    return 0;
    
exception_:
    perror("Program Failed Initialization");
    exit (EXIT_FAILURE);
}

/*
 Randomly generates a solution for simulated annealing or 
 foolish hill climbing and computes its fitness. 
 
 @param dst     A pointer to the structure that will contain 
                the randomly generated solution and its fitness. 
 */
void sima_rand (roulette_s *dst)
{
    uint32_t i;
    
    for (i = 0; i < pool_->solusize; i++) {
        ((uint16_t *)&dst->ptr[i])[0] = (uint16_t)rand();
        ((uint16_t *)&dst->ptr[i])[1] = (uint16_t)rand();
        ((uint16_t *)&dst->ptr[i])[2] = (uint16_t)rand();
        ((uint16_t *)&dst->ptr[i])[3] = (uint16_t)rand();
    }
    dst->ptr[i-1] &= pool_->cmask;
    dst->fitness = getfitness(dst->ptr);
}

/*
 The right side of the 'OR' statement in the simulated annealing algorithm: 
    if( h(NewS) < h(S) or random < e ^ [ (h(S) - h(NewS)) / T ] )
 This function is dynamically bound to pool_->e_pow. 
 
 @return    Returns the result of the boolean statement for the simulated annealing. 
 */
int e_pow_sa (void)
{
    return (SIMA_RAND() < pow (M_E, (pool_->rbuf[SIMA_curr].fitness - pool_->rbuf[SIMA_tmp].fitness) / pool_->T));
}

/*
 This replaces the the right side of the 'OR' statement in simulated annealing 
 with 'false', effectively removing it. This is for foolish hill climbing. 
 
 @return    Always returns 0, or 'false'. 
 */
int nop_hc (void)
{
    return 0;
}

/*
 Signal handler for SIGUSR1 (for child process). This simply determines 
 if simulated annealing/foolish hill climbing or the genetic algorithm
 is running, and calls the appropriate comand parser for whichever is running. 
 This signal is generated by the parent's command line, and sent to the child. 
 This parses the command buffer containing data piped to it from the parent. 
 
 @param signal  The integer value for the signal that invoked this handle.
 */
void cSIGUSR1 (int signal)
{
    int tmpint;
    unsigned char cbuf[CBUF_SIZE];
    gtoken_s *head;
    
    read(pipe_[0], cbuf, CBUF_SIZE);
    head = lex (cbuf);
    if (!head)
        printf ("Illegal Symbols Used\n");
    else {
        if (pool_->is_not_ge)
            csaparse();
        else
            cgeparse ();
        freetokens (stream_);
    }
    printf("> ");
    fflush (stdout);
    kill(getppid(), SIGUSR1);
}

/*
 Signal handler for SIGINT (for parent process). This function is 
 called when the user presses cntrl+c, and this sends the 'status'
 command to the child process by generating SIGUSR1 and piping 
 "status" to its command buffer. 
 
 @param signal  The integer value for the signal that invoked this handle.
 */
void pSIGINT (int signal)
{
    unsigned char cbuf[CBUF_SIZE];
    
    strcpy(cbuf, "status");
    cbuf[6] = UEOF;
    write(pipe_[1], cbuf, CBUF_SIZE);
    kill (pid_, SIGUSR1);
    pause();
}

/* 
 A NOP signal handler. This is used when SIGUSR1 is sent to parent
 process to wake it up from a pause. 
 */
void sigNOP (int signal){}

/***************** Functions To Print Output *****************/

void printgestatus (void)
{
    printf("Status at Generation: %llu\n", pool_->gen);
    printpool();
    printf("\nMost Fit ( weight = %f ):\n", pool_->rbuf[POOLSIZE-1].fitness);
    printchrom (pool_->rbuf[POOLSIZE-1].ptr);
    if (isfeasible(pool_->rbuf[POOLSIZE-1].ptr))
        printf("\nIs Feasible\n");
    else {
        printf("\nNot Feasible\n");
        if (pool_->bestfeasible) {
            printf("\nMost Fit Feasible ( weight = %f ):\n", getfitness(pool_->bestfeasible));
            printchrom (pool_->bestfeasible);
            printf("\n");
        }
        else
            printf("No Feasibles Found\n");
    }
    printf("Elapsed Time: %llu\n", (uint64_t)(time(NULL) - pool_->start));
    
}

void printsastatus (void)
{
    printf("Status at \"Generation\": %llu\n", pool_->gen);
    printchrom(pool_->rbuf[SIMA_curr].ptr);
    printf("\tFitness: %f\n", getfitness(pool_->rbuf[SIMA_curr].ptr));
    if (isfeasible(pool_->rbuf[SIMA_curr].ptr))
        printf ("Is Feasible\n");
    else {
        printf ("Is not Feasible\n");
        if (isfeasible(pool_->bestfeasible)) {
            printf("Best Feasible is: \n");
            printchrom(pool_->bestfeasible);
            printf("\tFitness: %f\n", getfitness(pool_->bestfeasible));
        }
        else
            printf ("No 'Good' Feasibles Found.\n");
    }
    printf("Elapsed Time: %llu\n", (uint64_t)(time(NULL) - pool_->start));
}


void printqword (uint64_t lword, uint8_t end)
{
    uint8_t i;
    
    if (!end)
        end = 64;
    for (i = 0; i < end; i++)
        printf("%llu",(lword >> i) & 1);
}

void printchrom (uint64_t *chrom)
{
    uint16_t i, n;
    
    n = pool_->chromsize;
    for (i = 0; i < n; i++) {
        if (i == n-1)
            printqword (chrom[i], pool_->remain);
        else
            printqword (chrom[i], 64);
    }
}

void printpool (void)
{
    uint16_t i, n;
    uint64_t *ptr;
    
    n = pool_->chromsize;
    ptr = pool_->popul;
    for (i = 0; i < POOLSIZE; i++) {
        printf("%d:\t", i);
        printchrom (ptr);
        printf("  %f, %d, %d", getfitness (ptr), ((uint8_t *)ptr)[7] >> (8 - pool_->bitlen), isfeasible (ptr));
        printf("\n");
        ptr += n;
    }
}

void insert_solset (solset_s **sset, vertex_s *v)
{
    solset_s *tmp;
    
    tmp = malloc (sizeof(*tmp));
    if (!tmp) {
        perror("Heap Allocation Error");
        exit(EXIT_FAILURE);
    }
    tmp->v = v;
    if (!*sset) {
        tmp->next = NULL;
        *sset = tmp;
    }
    else {
        tmp->next = *sset;
        *sset = tmp;
    }
}

void print_solset (solset_s *solset)
{
    int i;
    solset_s *iter, *tmp;
    
    printf("\n{");
    for (iter = solset, tmp = iter, i = 0; iter; tmp = iter, i++) {
        if (iter->next) {
            if (!(i % 8) && iter->next)
                printf("\n\t");
            printf ("%s,\t", iter->v->name);
        }
        else {
            if (!(i % 8))
                printf("\n\t");
            printf ("%s\t", iter->v->name);
        }
        iter = iter->next;
        free(tmp);
    }
    printf("\n}\n");
}

void printsolution (int index, uint64_t *ptr)
{
    int       i;
    float     fitness;
    uint64_t  *chrom;
    int       v1size, v2size;
    solset_s  *v1,    *v2;
    
    if (ptr)
        chrom = ptr;
    else
        chrom = pool_->rbuf[index].ptr;
    fitness = getfitness(chrom);
    if (isfeasible(chrom))
        printf ("Showing Feasible Chromosome %d with fitness %f:\n", index, fitness);
    else
        printf ("Showing Infeasible Chromosome %d with fitness %f:\n", index, fitness);
    printchrom(chrom);
    printf("\n");
    for (i = 0, v1size = 0, v2size = 0, v1 = NULL, v2 = NULL; i < pool_->bitlen; i++) {
        if (chrom[i / 64] & (1llu << (i % 64))) {
            ++v1size;
            insert_solset (&v1, pool_->graph->vtable[i]);
        }
        else {
            ++v2size;
            insert_solset (&v2, pool_->graph->vtable[i]);
        }
    }
    printf ("\nV1 = ");
    print_solset (v1);
    printf ("Length = %d\n\n", v1size);
    printf ("V2 = ");
    print_solset (v2);
    printf ("Length = %d\n", v2size);
}