/*
* Main GA Engine
* Author: Jonathan Hamm
*/

#include "gparse.h"
#include "bis.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <unistd.h>
#include <signal.h>

/*
 Counts the set bits in a long word:
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

/* Signal Handlers */
static void sigNOP (int signal) {}
static void cSIGUSR1 (int signal);
static void pSIGINT (int signal);
static void pSIGFPE (int signal);

static void pool_s_ (wgraph_s *g);
static void printqword (uint64_t lword, uint8_t mask);
static double sumweights (uint64_t *chrom);
static double getfitness (uint64_t *chrom);
static void computeprob (void);
static int isfeasible (uint64_t *chrom);
static void printchrom (uint64_t *chrom);

static pid_t pid_;
static int pipe_[2];
pool_s *pool_;

int prcmp (roulette_s *a, roulette_s *b)
{
  if (a->prob < b->prob)
    return -1;
  if (a->prob > b->prob)
    return 1;
  return 0;
}

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

int countdigits (uint64_t *cptr)
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

void printweights (void)
{
  uint64_t *ptr;
  uint16_t n, i;
  
  n = pool_->chromsize;
  ptr = pool_->popul;
  for (i = 0; i < POOLSIZE; i++) {
    printf("Weight: %f, %d, %d\n",sumweights (ptr), n, isfeasible(ptr));
    ptr += n;
  }
}

double getfitness (uint64_t *chrom)
{
  int setcount, differ;
  
  setcount = countdigits (chrom);
  differ = abs(2*setcount - pool_->bitlen);
  return sumweights (chrom) + (differ << 4);
}

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

int isfeasible (uint64_t *chrom)
{
  if (pool_->bitlen % 2)
    return  (2 * countdigits (chrom) == pool_->bitlen - 1)
            ||
            (2 * countdigits (chrom) == pool_->bitlen + 1);
  else
    return  (2 * countdigits (chrom) == pool_->bitlen);
}

uint8_t getbit (uint64_t *chrom, uint16_t pos)
{
  return (chrom[pos / 64] >> (pos % 64)) & 1llu;
}

void setbit (uint64_t *chrom, uint16_t pos, uint8_t val)
{
  chrom[pos / 64] &= ~(1llu << (pos % 64));
  if (val)
    chrom[pos / 64] |= (1llu << (pos % 64));
}

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
 that 'key' fits into.
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
 p1 p2  mask  | f
 ----------------
 0  0   0     | 0
 0  0   1     | 0
 0  1   0     | 0
 0  1   1     | 1
 1  0   0     | 1
 1  0   1     | 0
 1  1   0     | 1
 1  1   1     | 1
 
 Karnough Map
    00 01 11 10
 ---------------
 0|  0  0  1  1
 1|  0  1  1  0
 
 f = (~mask & p1) | p2
 
 Uniform Crossover
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

void mutate1 (uint64_t *victim)
{
  int i, n;
  uint16_t index;
  
  n = rand() % (pool_->bitlen / 8);
  if (!isfeasible(victim))
    ++n;
  for (i = 0; i < n; i++) {
    index = rand() % pool_->bitlen;
    victim[index / 64] ^= (1llu << (index % 64));
  }
}

void mutate2 (uint64_t *victim)
{
  uint16_t index, i;
  
  for (i = 0; i < 10; i++) {
    index = rand() % pool_->bitlen;
    victim[index / 64] ^= (1llu << (index % 64));
  }
}

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
    THROW_EXCEPTION();
  if (signal(SIGALRM, sigNOP) == SIG_ERR)
    THROW_EXCEPTION();
  memset (cbuf, 0, sizeof(cbuf));
  pool_s_ (g);
  if (!pool_)
    THROW_EXCEPTION();
  if (pipe (pipe_) == -1)
    THROW_EXCEPTION();
  pid_ = fork();
  if (pid_) {
    printf ("Now running program with chromosome size: %d\n> ", pool_->bitlen);
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
      pause ();
    }
  }
  else {
    if (signal(SIGUSR1, cSIGUSR1) == SIG_ERR)
      THROW_EXCEPTION();
    if (signal(SIGINT, sigNOP) == SIG_ERR)
      THROW_EXCEPTION();
    signal(SIGQUIT, pSIGFPE);
    n = pool_->chromsize;
    pool_->start = clock();
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
          for (rpt = &pool_->rbuf[rand()%POOLSIZE];
               rpt->ptr == pool_->bestfeasible ||
               rpt->ptr == dst1 ||
               rpt->ptr == dst2
               ;rpt = &pool_->rbuf[rand()%POOLSIZE]);
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
      ++pool_->gen;
    }
  }
  return 0;
  
exception_:
  perror("Failed Setting Program Initialization");
  exit (EXIT_FAILURE);
}

void printstatus (void)
{
  int i;
  
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
    }
    else
      printf("No Feasibles Found\n");
  }
  printf("Elapsed Time: %llu\n", (uint64_t)(clock() - pool_->start));

}


void cSIGUSR1 (int signal)
{
  int tmpint;
  unsigned char cbuf[CBUF_SIZE];
  gtoken_s *head;
  
  read(pipe_[0], cbuf, CBUF_SIZE);
  head = lex_ (cbuf);
  if (!head)
    printf ("Illegal Symbols Used\n");
  else {
    cparse ();
    freetokens (stream_);
  }
  printf("> ");
  fflush (stdout);
  kill(getppid(), SIGALRM);
}

void pSIGINT (int signal)
{
  unsigned char cbuf[CBUF_SIZE];
  
  strcpy(cbuf, "status");
  cbuf[6] = UEOF;
  write(pipe_[1], cbuf, CBUF_SIZE);
  kill (pid_, SIGUSR1);
  pause();
}

void pSIGFPE (int signal)
{
  kill (getppid(), SIGSEGV);
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

void printsolution (int index)
{
  int       i;
  float     fitness;
  uint64_t  *chrom;
  int       v1size, v2size;
  solset_s  *v1,    *v2;
  
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