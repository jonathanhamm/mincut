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

/* Signal Handlers */
static void sigNOP (int signal) {}
static void cSIGUSR1 (int signal);
static void pSIGINT (int signal);

static pool_s *pool_s_ (wgraph_s *g);
static void printqword (uint64_t lword, uint8_t mask);
static float sumweights (pool_s *p, uint64_t *chrom);
static float getfitness (pool_s *p, uint64_t *chrom);
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
  float sum;
  
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
  p->select = roulette_sf;
  p->cross = singlepoint_cr;
  p->mutate = mutate1;
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
    p->accum += p->rbuf[i].prob;
    p->rbuf[i].cummulative = (int)p->rbuf[i].prob;
    if (i)
      p->rbuf[i].cummulative += p->rbuf[i-1].cummulative;
  }
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

/*
 Iterative Binary search that returns the index of roulette 'region'
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

void roulette_sf (pool_s *p, selected_s *parents)
{
  int i;
  uint32_t i1, i2;
  
  for (i = 0; i < NSELECT; i++) {
    i1 = bsearch_r (p->rbuf, rand() % p->accum);
    i2 = bsearch_r (p->rbuf, rand() % p->accum);
    if (i2 == i1)
      i2 = (i2 - 1) % POOLSIZE;
    parents->couples[i].p1 = &p->rbuf[i1];
    parents->couples[i].p2 = &p->rbuf[i2];
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

float sumweights (pool_s *p, uint64_t *chrom)
{
  uint8_t pos;
  uint16_t i, j, csize, nedges;
  uint64_t iter,
           *ptr;
  vertex_s *v;
  edge_s **edges;
  float weight;

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

float getfitness (pool_s *p, uint64_t *chrom)
{
  int setcount, differ;
  
  setcount = countdigits (p, chrom);
  differ = abs(2*setcount-GET_CHBITLEN(p));
  return sumweights (p, chrom) + (differ << 4);
}

void computeprob (pool_s *p)
{
  uint16_t i, n;
  uint64_t *ptr;
  float sum, fc1, fc2;
  
  sum = p->fitsum;
  for (i = 0, p->accum = 0; i < POOLSIZE; i++)
    p->rbuf[i].prob = sum / p->rbuf[i].fitness;
  qsort (p->rbuf, POOLSIZE, sizeof(roulette_s), (int (*)(const void *, const void *))prcmp);
  for (i = 0; i < POOLSIZE; i++) {
    p->rbuf[i].cummulative = (int)p->rbuf[i].prob;
    if (i)
      p->rbuf[i].cummulative += p->rbuf[i-1].cummulative;
    p->accum += p->rbuf[i].prob;
  }
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

int getbit (uint64_t *chrom, uint16_t pos)
{
  return (chrom[pos / 64] >> (pos % 64)) & 1llu;
}

void setbit (uint64_t *chrom, uint16_t pos, int val)
{
  chrom[pos / 64] &= ~(1llu << (pos % 64));
  if (val)
    chrom[pos / 64] |= (1llu << (pos % 64));
}

#define TESTCROSSOVER_ 
#undef TESTCROSSOVER_

void singlepoint_cr (pool_s *p, roulette_s *rp1, roulette_s *rp2)
{
  int i, pindex;
  uint16_t point, inv;
  uint64_t mask;
  uint64_t *backup1, *backup2, *p1, *p2,
           *dst1, *dst2;
  
  p1 = rp1->ptr;
  p2 = rp2->ptr;
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
    rp2 = &p->rbuf[1%POOLSIZE];
    dst2 = rp2->ptr;
  }
  p->fitsum -= (rp1->fitness + rp2->fitness);
#ifdef TESTCROSSOVER_
  printf ("parents: \n");
  printchrom(p, p1);
  printf("\n");
  printchrom (p, p2);
#endif
  backup1 = &p->popul[CRBACKUP1];
  backup2 = &p->popul[CRBACKUP2];
  for (i = 0; i < p->chromsize; i++)
    backup1[i] = p1[i];
  for (i = 0; i < p->chromsize; i++)
    backup2[i] = p2[i];
  point = rand() % GET_CHBITLEN(p);
  inv = GET_CHBITLEN(p) - point;
#ifdef TESTCROSSOVER_
  printf("\nPoint: %d\n", point);
#endif
  if (point >= inv) {
    for (i = 0; i < point; i++) {
      setbit(dst1, i, getbit(backup2, point+i));
      setbit(dst2, inv+i, getbit(backup1, i));
    }
  }
  else {
    for (i = 0; i < point; i++) {
      setbit(dst1, i, getbit(backup2, inv+i));
      setbit(dst2, point+i, getbit(backup1, i));
    }

  }
#ifdef TESTCROSSOVER_
  printf("Children:\n");
  printchrom(p, dst1);
  printf("\n");
  printchrom(p, dst2);
  printf("\n\n");
#endif
  rp1->fitness = getfitness (p, dst1);
  rp2->fitness = getfitness (p, dst2);
  p->fitsum += (rp1->fitness + rp2->fitness);
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
void uniform_cr (pool_s *p, roulette_s *rp1, roulette_s *rp2)
{
  uint16_t i;
  uint64_t  mask,  backup,
            *dst1,  *dst2,
            *p1,    *p2;
  
  p1 = rp1->ptr;
  p2 = rp2->ptr;
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
    rp2 = &p->rbuf[1%POOLSIZE];
    dst2 = rp2->ptr;
  }
  p->fitsum -= (rp1->fitness + rp2->fitness);
  for (i = 0; i < p->chromsize; i++) {
    ((uint16_t *)&mask)[0] = (uint16_t)rand();
    ((uint16_t *)&mask)[1] = (uint16_t)rand();
    ((uint16_t *)&mask)[2] = (uint16_t)rand();
    ((uint16_t *)&mask)[3] = (uint16_t)rand();
    backup = p1[i];
    dst1[i] = (~mask & backup) | (mask & p2[i]);
    dst2[i] = (mask & backup) | (~mask & p2[i]);
  }
  if (rand() % 100 < p->mutateprob) {
    p->mutate (p, dst1);
    p->mutate (p, dst2);
  }
  rp1->fitness = getfitness (p, dst1);
  rp2->fitness = getfitness (p, dst2);
  p->fitsum += (rp1->fitness + rp2->fitness);
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
    n = p->chromsize;
    p->start = clock();
    while (1) {
      p->select (p, &parents);
      for (i = 0; i < NSELECT; i++)
        p->cross (p, parents.couples[i].p1, parents.couples[i].p2);
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

#define COM_PROB 1
#define COM_OP   2
static void cparse (gtoken_s *stream);
static int  p_mutate (gtoken_s *stream);
static int  p_mutate_ (gtoken_s *stream);
static int  p_mval (gtoken_s *stream);
static void p_show (gtoken_s *stream);
static void p_feasible (gtoken_s *stream);


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
    cparse (head);
  freetokens (head);
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

void printsolution (pool_s *p, int index)
{
  
}

inline void printfittest (pool_s *p)
{
  printsolution (p, POOLSIZE-1);
}

void cparse (gtoken_s *stream)
{
  int result, val;
  
  if (
      !strcmp (stream->lexeme, "exit")  ||
      !strcmp (stream->lexeme, "quit")  ||
      !strcmp (stream->lexeme, "Quit")  ||
      !strcmp (stream->lexeme, "q")     ||
      !strcmp (stream->lexeme, "Q")
      )
  {
    printf("Final:\n");
    printstatus ();
    kill(getppid(), SIGQUIT);
    exit(EXIT_SUCCESS);
  }
  else if (!strcmp(stream->lexeme, "set")) {
    result = p_mutate (stream->next);
    val = p_mval (stream->next);
    if (result == COM_PROB) {
      if (val < 0)
        val = 0;
      else if (val > 100)
        val = 100;
      pool_->mutateprob = (uint8_t)val;
      printf ("Mutation Probability Set to: %d\n", val);
    }
    else if (result == COM_OP) {
      if (val <= 1) {
        pool_->mutate = mutate1;
        printf ("Using Mutation Operator 1\n");
      }
      else {
        pool_->mutate = mutate2;
        printf ("Using Mutation Operator 2\n");
      }
    }
  }
  else if (!strcmp(stream->lexeme, "get")) {
    result = p_mutate (stream->next);
    if (result == COM_PROB)
      printf ("Current Mutation Probability: %d\n", pool_->mutateprob);
    else if (result == COM_OP)
      printf ("Currently Using Mutation Operator: %d\n", (pool_->mutate == mutate1) ? 1 : 2);
  }
  else if (!strcmp(stream->lexeme, "show")) {
    p_show (stream->next);
  }
  else if (!strcmp (stream->lexeme, "status")) {
    printstatus ();
  }
  else {
    printf ("Command Line Error: Unrecognized: '%s'\n", stream->lexeme);
  }
}

int p_mutate (gtoken_s *stream)
{
  if (!strcmp (stream->lexeme, "mutate"))
    return p_mutate_ (stream->next);
  else
    printf ("Command Line Error: Unrecognized: '%s'\n", stream->lexeme);
}

int p_mutate_ (gtoken_s *stream)
{
  if (!strcmp (stream->lexeme, "prob"))
    return COM_PROB;
  else if (!strcmp (stream->lexeme, "op"))
    return COM_OP;
  printf ("Command Line Error: Unrecognized: '%s'\n", stream->lexeme);
  return 0;
}

int p_mval (gtoken_s *stream)
{
  if (stream->type == T_NUM)
    return atoi (stream->lexeme);
  else
    printf ("Expected Number, but got %s\n", stream->lexeme);
}

void p_show (gtoken_s *stream)
{
  if (stream->type == T_NUM) {
    
  }
  else if (!strcmp (stream->lexeme, "best")) {
    p_feasible (stream->next);
  }
  else
    printf ("Expected number or 'best', but got %s\n", stream->lexeme);

}

void p_feasible (gtoken_s *stream)
{
  if (!strcmp (stream->lexeme, "feasible")) {
    
  }
  else if (stream->type != T_EOF) {
    printf ("Expected 'feasible' or nothing, but got %s\n", stream->lexeme);
  }
}