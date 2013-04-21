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

static pool_s *pool_s_ (wgraph_s *g);
static void printqword (uint64_t lword, uint8_t mask);
static double sumweights (pool_s *p, uint64_t *chrom);
static double getfitness (pool_s *p, uint64_t *chrom);
static void computeprob (pool_s *p);
static int isfeasible (pool_s *p, uint64_t *chrom);
static void printchrom (pool_s *p, uint64_t *chrom);

static pid_t pid_;
static int pipe_[2];
static pool_s *pool_;

int prcmp (roulette_s *a, roulette_s *b)
{
  if (a->prob < b->prob)
    return -1;
  if (a->prob > b->prob)
    return 1;
  return 0;
}

pool_s *pool_s_ (wgraph_s *g)
{
  pool_s *p;
  uint16_t  i, j;
  uint64_t  *ptr;
  double sum;
  
  p = calloc(1, sizeof(*p) + POOLSIZE * CQWORDSIZE(g->nvert+1) * 8);
  if (!p) {
    perror ("Malloc Error");
    exit (EXIT_FAILURE);
  }
  pool_ = p;
  p->chromsize = (g->nvert / 64) + (g->nvert % 64 != 0);
  p->remain = g->nvert % 64;
  if (!(g->nvert % 64) && g->nvert)
    p->cmask = 0xffffffffffffffffllu;
  else {
    for (i = 0; i < p->remain; i++)
        p->cmask |= (1llu << i);
  }
  p->select = tournament_sf;
  p->cross = uniform_cr;
  p->mutate = mutate1;
  p->k = TOURN_K;
  p->gen = 0;
  ptr = p->popul;
  for (i = 0; i < POOLSIZE; i++) {
    for (j = 0; j < p->chromsize; j++, ptr++) {
      ((uint16_t *)ptr)[0] = (uint16_t)rand();
      ((uint16_t *)ptr)[1] = (uint16_t)rand();
      ((uint16_t *)ptr)[2] = (uint16_t)rand();
      ((uint16_t *)ptr)[3] = (uint16_t)rand();
    }
    *(ptr-1) &= p->cmask;
  }
  p->graph = g;
  p->mutateprob = INITMUTATIONPROB;
  p->bitlen = ((p->remain) ? ((p->chromsize * 64) - (64 - p->remain)) : p->chromsize*64);
  for (sum = 0, ptr = p->popul, i = 0; i < POOLSIZE; i++, ptr += p->chromsize) {
    p->rbuf[i].fitness = getfitness (p, ptr);
    p->rbuf[i].ptr = ptr;
    sum += p->rbuf[i].fitness;
  }
  for (i = 0, p->accum = 0; i < POOLSIZE; i++)
    p->rbuf[i].prob = sum / p->rbuf[i].fitness;
  qsort (p->rbuf, POOLSIZE, sizeof(roulette_s), (int (*)(const void *, const void *))prcmp);
  for (i = 0; i < POOLSIZE; i++) {
    p->accum += (int)p->rbuf[i].prob;
    p->rbuf[i].cummulative = (i+1);
    if (i)
      p->rbuf[i].cummulative += p->rbuf[i-1].cummulative;
  }
  p->ranksum = p->rbuf[POOLSIZE-1].cummulative;
  p->fitsum = sum;
  return p;
}

void printqword (uint64_t lword, uint8_t end)
{
  uint8_t i;
  
  if (!end)
    end = 64;
  for (i = 0; i < end; i++)
    printf("%llu",(lword >> i) & 1);
}

void printchrom (pool_s *p, uint64_t *chrom)
{
  uint16_t i, n;
  
  n = p->chromsize;
  for (i = 0; i < n; i++) {
    if (i == n-1)
      printqword (chrom[i], p->remain);
    else
      printqword (chrom[i], 64);
  }
}

void printpool (pool_s *p)
{
  uint16_t i, n, index;
  uint64_t *ptr;
  
  n = p->chromsize;
  ptr = p->popul;
  for (index = 0; index < POOLSIZE; index++) {
    printf("%d:\t", index);
    printchrom (p, ptr);
    printf("  %f, %d, %d", getfitness (p, ptr), ((uint8_t *)ptr)[7] >> (8 - GET_CHBITLEN(p)), isfeasible (p, ptr));
    printf("\n");
    ptr += n;
  }
}

int countdigits (pool_s *p, uint64_t *cptr)
{
  int i, count, n;

  n = p->chromsize-1;
  for (i = 0, count = 0; i < n; i++) {
    count += countbitsLW((uint32_t)cptr[i]);
    count += countbitsLW((uint32_t)(cptr[i] >> 32));
  }
  count += countbitsLW((uint32_t)(cptr[i] & p->cmask));
  count += countbitsLW((uint32_t)((cptr[i] & p->cmask) >> 32));
  return count;
}

int iscut(pool_s *p, uint64_t *chrom, vertex_s *v)
{
  uint16_t i, j, nvert;
  uint64_t iter;
  vertex_s **vtable;
  
  nvert = p->graph->nvert;
  vtable = p->graph->vtable;
  i = vgetindex (v);
  if ((~chrom[i / 64] & (1llu << (i % 64))))
    return 1;
  return 0;
}

double sumweights (pool_s *p, uint64_t *chrom)
{
  uint8_t pos;
  uint16_t i, j, csize, nedges;
  uint64_t iter,
           *ptr;
  vertex_s *v;
  edge_s **edges;
  double weight;

  for (weight = 0, i = 0, ptr = chrom; i < p->chromsize; i++, ptr++) {
    if (i == p->chromsize-1 && p->remain)
      csize = p->remain;
    else
      csize = 64;
    for (iter = *ptr, pos = 0; pos < csize; iter &= ~(1llu << pos)) {
      pos = ffsl(iter);
      if (!pos)
        break;
      --pos;
      v = p->graph->vtable[/* i*64 */ (i << 6) + pos];
      nedges = v->nedges;
      edges = v->edges;
      for (j = 0; j < nedges; j++) {
        if (iscut(p, chrom, (edges[j]->v1 == v) ? edges[j]->v2 : edges[j]->v1))
          weight += edges[j]->weight;
      }
    }
  }
  return weight;
}

void printweights (pool_s *p)
{
  uint64_t *ptr;
  uint16_t n, i;
  
  n = p->chromsize;
  ptr = p->popul;
  for (i = 0; i < POOLSIZE; i++) {
    printf("Weight: %f, %d, %d\n",sumweights (p, ptr), n, isfeasible(p, ptr));
    ptr += n;
  }
}

double getfitness (pool_s *p, uint64_t *chrom)
{
  int setcount, differ;
  
  setcount = countdigits (p, chrom);
  differ = abs(2*setcount - GET_CHBITLEN(p));
  return sumweights (p, chrom) + (differ << 4);
}

void computeprob (pool_s *p)
{
  uint16_t i, n;
  double sum;
  
  sum = p->fitsum;
  for (p->accum = 0, i = 0; i < POOLSIZE; i++) {
    p->rbuf[i].prob = sum / p->rbuf[i].fitness;
    p->accum += (int)p->rbuf[i].prob;
  }
  qsort (p->rbuf, POOLSIZE, sizeof(roulette_s), (int (*)(const void *, const void *))prcmp);
  if (isfeasible(p, p->rbuf[POOLSIZE-1].ptr))
    p->bestfeasible = p->rbuf[POOLSIZE-1].ptr;
}

int isfeasible (pool_s *p, uint64_t *chrom)
{
  if (GET_CHBITLEN(p) % 2)
    return  (2 * countdigits (p, chrom) == GET_CHBITLEN(p) - 1)
            ||
            (2 * countdigits (p, chrom) == GET_CHBITLEN(p) + 1);
  else
    return  (2 * countdigits (p, chrom) == GET_CHBITLEN(p));
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

void roulette_sf (pool_s *p, selected_s *parents)
{
  int i, j;
  uint32_t i1, i2;
  float sum;
  
  for (i = 0; i < NSELECT; i++) {
    i1 = rand() % p->accum;
    i2 = rand() % p->accum;
    for (j = 0, sum = 0; j < POOLSIZE && sum <= (float)i1; sum += p->rbuf[j].prob, j++);
    i1 = (j == POOLSIZE) ? j-1 : j;
    for (j = 0, sum = 0; j < POOLSIZE && sum <= (float)i2; sum += p->rbuf[j].prob, j++);
    i2 = (j == POOLSIZE) ? j-1 : j;
    if (i2 == i1)
      i2 = (i2 - 1) % POOLSIZE;
    parents->couples[i].p1 = &p->rbuf[i1];
    parents->couples[i].p2 = &p->rbuf[i2];
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

void rank_sf (pool_s *p, selected_s *parents)
{
  int i;
  uint32_t i1, i2;
  for (i = 0; i < NSELECT; i++) {
    i1 = bsearch_r (p->rbuf, rand() % p->ranksum);
    i2 = bsearch_r (p->rbuf, rand() % p->ranksum);
    if (i2 == i1)
      i2 = (i2 - 1) % POOLSIZE;
    parents->couples[i].p1 = &p->rbuf[i1];
    parents->couples[i].p2 = &p->rbuf[i2];
  }
}

void tournament_sf (pool_s *p, selected_s *parents)
{
  int       i;
  uint8_t   k, R;
  uint32_t  i1a, i1b,
            i2a, i2b;
  
  for (i = 0, k = p->k; i < NSELECT; i++) {
    i1a = rand() % POOLSIZE;
    while ((i1b = rand() % POOLSIZE) == i1a);
    i2a = rand() % POOLSIZE;
    while ((i2b = rand() % POOLSIZE) == i2a);
    R  = rand() % 100;
    if (R < k) {
      if (i1a > i1b)
        parents->couples[i].p1 = &p->rbuf[i1a];
      else
        parents->couples[i].p1 = &p->rbuf[i1b];
      if (i2a > i2b)
        parents->couples[i].p2 = &p->rbuf[i2a];
      else
        parents->couples[i].p2 = &p->rbuf[i2b];
    }
    else {
      if (i1a < i1b)
        parents->couples[i].p1 = &p->rbuf[i1a];
      else
        parents->couples[i].p1 = &p->rbuf[i1b];
      if (i2a < i2b)
        parents->couples[i].p2 = &p->rbuf[i2a];
      else
        parents->couples[i].p2 = &p->rbuf[i2b];
    }
  }
}

void npoint_cr (pool_s *p, uint64_t *p1, uint64_t *p2, uint64_t *dst1, uint64_t *dst2)
{
  int i, j, pindex;
  uint16_t point, inv;
  uint64_t mask;
  uint64_t *backup1, *backup2;
  
  backup1 = &p->popul[CRBACKUP1];
  backup2 = &p->popul[CRBACKUP2];
  for (i = 0; i < p->chromsize; i++)
    backup1[i] = p1[i];
  for (i = 0; i < p->chromsize; i++)
    backup2[i] = p2[i];
  point = rand() % (GET_CHBITLEN(p) / CR_N);
  inv = GET_CHBITLEN(p) - point;
  for (i = 0; i < CR_N; i++) {
    for (j = i * (GET_CHBITLEN(p) / CR_N); j < inv; j++) {
      setbit(dst1, j, getbit(backup2, point+j));
      setbit(dst2, j, getbit(backup1, point+j));
    }
    for (j = i * (GET_CHBITLEN(p) / CR_N); j < point; j++) {
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
void uniform_cr (pool_s *p, uint64_t *p1, uint64_t *p2, uint64_t *dst1, uint64_t *dst2)
{
  uint16_t i;
  uint64_t  mask,  backup;

  for (i = 0; i < p->chromsize; i++) {
    ((uint16_t *)&mask)[0] = (uint16_t)rand();
    ((uint16_t *)&mask)[1] = (uint16_t)rand();
    ((uint16_t *)&mask)[2] = (uint16_t)rand();
    ((uint16_t *)&mask)[3] = (uint16_t)rand();
    backup = p1[i];
    dst1[i] = (~mask & backup) | (mask & p2[i]);
    dst2[i] = (mask & backup) | (~mask & p2[i]);
  }
}

void mutate1 (pool_s *p, uint64_t *victim)
{
  int i, n;
  uint16_t index;
  
  n = rand() % (GET_CHBITLEN(p) / 8);
  if (!isfeasible(p, victim))
    ++n;
  for (i = 0; i < n; i++) {
    index = rand() % GET_CHBITLEN(p);
    victim[index / 64] ^= (1llu << (index % 64));
  }
}

void mutate2 (pool_s *p, uint64_t *victim)
{
  uint16_t index, i;
  
  for (i = 0; i < 10; i++) {
    index = rand() % GET_CHBITLEN(p);
    victim[index / 64] ^= (1llu << (index % 64));
  }
}

int run_ge (wgraph_s *g)
{
  pool_s *p;
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
  p = pool_s_ (g);
  if (!p)
    THROW_EXCEPTION();
  if (pipe (pipe_) == -1)
    THROW_EXCEPTION();
  pid_ = fork();
  if (pid_) {
    printf ("Now running program.\n> ");
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
    n = p->chromsize;
    p->start = clock();
    while (1) {
      p->select (p, &parents);
      for (i = 0; i < NSELECT; i++) {
        rp1 = parents.couples[i].p1;
        rp2 = parents.couples[i].p2;
        p1  = rp1->ptr;
        p2  = rp2->ptr;
        if (rp1->ptr == p->bestfeasible) {
          rp1 = &p->rbuf[0];
          dst1 = rp1->ptr;
        }
        else
          dst1 = p1;
        if (rp2->ptr == p->bestfeasible) {
          rp2 = &p->rbuf[0];
          dst2 = rp2->ptr;
        }
        else
          dst2 = p2;
        if (dst1 == dst2) {
          for (rpt = &p->rbuf[rand()%POOLSIZE];
               rpt->ptr == p->bestfeasible ||
               rpt->ptr == dst1 ||
               rpt->ptr == dst2
               ;rpt = &p->rbuf[rand()%POOLSIZE]);
          rp2 = rpt;
          dst2 = rpt->ptr;
        }
        p->fitsum -= (rp1->fitness + rp2->fitness);
        p->cross (p, p1, p2, dst1, dst2);
        if (rand() % 100 < p->mutateprob) {
          p->mutate (p, dst1);
          p->mutate (p, dst2);
        }
        rp1->fitness = getfitness (p, dst1);
        rp2->fitness = getfitness (p, dst2);
        p->fitsum += (rp1->fitness + rp2->fitness);

      }
      computeprob(p);
      ++p->gen;
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
  printpool(pool_);
  printf("\nMost Fit ( weight = %f ):\n", pool_->rbuf[POOLSIZE-1].fitness);
  printchrom (pool_, pool_->rbuf[POOLSIZE-1].ptr);
  if (isfeasible(pool_, pool_->rbuf[POOLSIZE-1].ptr))
    printf("\nIs Feasible\n");
  else {
    printf("\nNot Feasible\n");
    if (pool_->bestfeasible) {
      printf("\nMost Fit Feasible ( weight = %f ):\n", getfitness(pool_, pool_->bestfeasible));
      printchrom (pool_, pool_->bestfeasible);
    }
    else
      printf("No Feasibles Found\n");
  }
  printf("Elapsed Time: %llu\n", (uint64_t)(clock() - pool_->start));

}

#define COM_ERR  0
#define COM_PROB 1
#define COM_OP   2
#define COM_X    3
#define COM_SEL  4
static void cparse (void);
static int  p_op (void);
static int  p_mutate (void);
static void p_show (void);
static int p_feasible (void);

void cSIGUSR1 (int signal)
{
  int tmpint;
  unsigned char cbuf[CBUF_SIZE];
  gtoken_s *head;
  
  read(pipe_[0], cbuf, CBUF_SIZE);
  head = lex_ (cbuf);
  if (!head) {
    printf ("Illegal Symbols Used\n");
    goto exit_;
  }
  else
    cparse ();
  freetokens (stream_);
exit_:
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
    perror("Malloc Error\n");
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

void printsolution (pool_s *p, int index)
{
  int       i;
  float     fitness;
  uint64_t  *chrom;
  int       v1size, v2size;
  solset_s  *v1,    *v2,
            *iter,  *tmp;
  
  chrom = p->rbuf[index].ptr;
  fitness = getfitness(pool_, chrom);
  if (isfeasible(p, chrom))
    printf ("Showing Feasible Chromosome %d with fitness %f:\n", index, fitness);
  else 
    printf ("Showing Infeasible Chromosome %d with fitness %f:\n", index, fitness);
  printchrom(p, chrom);
  printf("\n");
  for (i = 0, v1size = 0, v2size = 0, v1 = NULL, v2 = NULL; i < GET_CHBITLEN(p); i++) {
    if (chrom[i / 64] & (1llu << (i % 64))) {
      ++v1size;
      insert_solset (&v1, p->graph->vtable[i]);
    }
    else {
      ++v2size;
      insert_solset (&v2, p->graph->vtable[i]);
    }
  }
  printf ("\nV1 = ");
  print_solset (v1);
  printf ("\nLength = %d\n", v1size);
  printf ("V2 = ");
  print_solset (v2);
  printf ("\nLength = %d\n", v2size);
}

void cparse (void)
{
  int result, val;
  
  if (
      !strcmp (stream_->lexeme, "exit")  ||
      !strcmp (stream_->lexeme, "quit")  ||
      !strcmp (stream_->lexeme, "Quit")  ||
      !strcmp (stream_->lexeme, "q")     ||
      !strcmp (stream_->lexeme, "Q")
      )
  {
    printf("Final:\n");
    printstatus ();
    kill(getppid(), SIGQUIT);
    exit(EXIT_SUCCESS);
  }
  else if (!strcmp(stream_->lexeme, "set")) {
    GTNEXT();
    result = p_op ();
    if (stream_->type == T_NUM) {
      val = atoi (stream_->lexeme);
      GTNEXT();
      switch (result) {
        case COM_OP:
          if (val <= 1) {
            pool_->mutate = mutate1;
            printf("Mutate operator now set to: 1\n");
          }
          else {
            pool_->mutate = mutate2;
            printf("Mutate operator now set to: 2\n");
          }
          break;
        case COM_PROB:
          if (val < 0) val = 0;
          else if (val > 100) val = 100;
          pool_->mutateprob = (uint8_t)val;
          printf("Mutate probability now set to %d\n", val);
          break;
        case COM_X:
          if (val <= 1) {
            pool_->cross = uniform_cr;
            printf("Now using uniform crossover\n");
          }
          else {
            pool_->cross = npoint_cr;
            printf("Now using n-point crossover, n = %d\n", CR_N);
          }
          break;
        case COM_SEL:
          if (val <= 1) {
            pool_->select = roulette_sf;
            printf("Now using roulette selection.");
          }
          else if (val == 2) {
            pool_->select = rank_sf;
            printf("Now using rank selection\n");
          }
          else {
            pool_->select = tournament_sf;
            printf("Now using tournament selection\n");
          }
          break;
        default:
          break;
      }
    }
    else
      printf ("Command Line Error: Expected number, but got '%s'\n", stream_->lexeme);
  }
  else if (!strcmp(stream_->lexeme, "get")) {
    GTNEXT();
    result = p_op();
    switch (result) {
      case COM_OP:
        if (pool_->mutate == mutate1)
          printf("Currently using mutation operator 1.\n");
        else
          printf("Currently using mutation operator 2.\n");
        break;
      case COM_PROB:
        printf("Current mutation probability is %u%%\n", pool_->mutateprob);
        break;
      case COM_X:
        if (pool_->cross == uniform_cr)
          printf("Currently using uniform crossover.\n");
        else
          printf ("Currently using n-point crossover, n = %d\n", CR_N);
        break;
      case COM_SEL:
        if (pool_->select == roulette_sf)
          printf("Currently using roulette selection.\n");
        else if (pool_->select == rank_sf)
          printf("Currently using rank selection.\n");
        else
          printf("Currently using tournament selection.\n");
        break;
      default:
        break;
    }

  }
  else if (!strcmp(stream_->lexeme, "show")) {
    GTNEXT();
    p_show();
  }
  else if (!strcmp (stream_->lexeme, "status")) {
    GTNEXT();
    printstatus();
  }
  else
    printf ("Command Line Error: Unrecognized: '%s'\n", stream_->lexeme);
}

int p_op (void)
{
  if (!strcmp(stream_->lexeme, "mutate")) {
    GTNEXT();
    return p_mutate();
  }
  else if (!strcmp(stream_->lexeme, "cross")) {
    GTNEXT();
    return COM_X;
  }
  else if (!strcmp(stream_->lexeme, "select")) {
    GTNEXT();
    return COM_SEL;
  }
  printf ("Command Line Error: Expected 'mutate', 'cross', or 'select', but got '%s'.\n", stream_->lexeme);
  return COM_ERR;
}

int p_mutate (void)
{
  if (!strcmp (stream_->lexeme, "op")) {
    GTNEXT();
    return COM_OP;
  }
  else if (!strcmp(stream_->lexeme, "prob")) {
    GTNEXT();
    return COM_PROB;
  }
  printf ("Command Line Error: Expected 'op' or 'prob' but got: '%s'.\n", stream_->lexeme);
  return COM_ERR;
}

void p_show (void)
{
  int result, index;
  
  if (stream_->type == T_NUM) {
    index = atoi(stream_->lexeme);
    GTNEXT();
    if (index <= 0 || index > POOLSIZE)
      printf ("Value %d out of range. Range is 1 to %d.\n", index, POOLSIZE);
    else
      printsolution(pool_, --index);
  }
  else if (!strcmp (stream_->lexeme, "best")) {
    GTNEXT();
    result = p_feasible ();
    if (!result)
      printsolution (pool_, POOLSIZE-1);
    else if (result == 1) {
      for (index = POOLSIZE-1; index >= 0
          && pool_->rbuf[index].ptr != pool_->bestfeasible; index--);
      printsolution (pool_, index);
    }

  }
  else
    printf ("Expected number or 'best', but got '%s'.\n", stream_->lexeme);
}

int p_feasible (void)
{
  int index;
  
  if (!strcmp (stream_->lexeme, "feasible")) {
    GTNEXT();
    return 1;
  }
  else if (stream_->type != T_EOF) {
    printf ("Expected 'feasible' or nothing, but got '%s'.\n", stream_->lexeme);
    return -1;
  }
  return 0;
}