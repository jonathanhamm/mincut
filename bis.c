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

static pool_s *pool_init (wgraph_s *g);

static void sigdummy (int signal);
static void siginthandle(int signal);
static void sigquithandle (int signal);

static pool_s *pool_s_ (uint16_t csize);
static void printlword (uint64_t lword, uint8_t mask);
static float sumweights (pool_s *p, uint64_t *chrom);
static float getfitness (pool_s *p, uint64_t *chrom);
static float computeprob (pool_s *p);
static int isfeasible (pool_s *p, uint64_t *chrom);
static void printchrom (pool_s *p, uint64_t *chrom);

static pool_s *pool_;

pool_s *pool_s_ (uint16_t csize)
{
  int i;
  pool_s *pool;
  
  pool = calloc(1, sizeof(*pool) + _POOLSIZE*_CQWORDSIZE(csize)*8);
  if (!pool)
    goto err_;
  pool_ = pool;
  pool->chromsize = (csize >> 6) + (csize % 64 != 0);
  pool->remain = csize % 64;
  for (i = 0; i < pool->remain; i++)
    pool->cmask |= (1 << i);
  pool->cross = singlepoint_cr;
  pool->mutate = mutate1;
  printf("QWORD size: %d\n", pool->chromsize);
  printf("Remainder: %d\nMask\n", pool->remain);
  printlword(pool->cmask, 64);
  printf ("\n");
  return pool;
  
err_:
  return NULL;
}

pool_s *pool_init (wgraph_s *g)
{
  uint16_t  i, j;
  pool_s    *pool;
  uint64_t  *ptr;
  
  pool = pool_s_(g->nvert);
  if (!pool)
    return NULL;
  pool->gen = 0;
  ptr = pool->popul;
  for (i = 0; i < _POOLSIZE; i++) {
    for (j = 0; j < pool->chromsize; j++, ptr++) {
      ((uint16_t *)ptr)[0] = (uint16_t)rand();
      ((uint16_t *)ptr)[1] = (uint16_t)rand();
      ((uint16_t *)ptr)[2] = (uint16_t)rand();
      ((uint16_t *)ptr)[3] = (uint16_t)rand();
      *ptr &= pool->cmask;
    }
  }
  pool->graph = g;
  return pool;
}

void printlword (uint64_t lword, uint8_t end)
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
      printlword (chrom[i], p->remain);
    else
      printlword (chrom[i], 64);
  }
}

void printpool (pool_s *p)
{
  uint16_t i, n, index;
  uint64_t *ptr;
  
  n = p->chromsize;
  ptr = p->popul;
  for (index = 0; index < _POOLSIZE; index++) {
    printf("%d:\t", index);
    printchrom (p, ptr);
    printf("  %f, %d, %d", getfitness (p, ptr), ((uint8_t *)ptr)[7] >> (8 - _GET_CHBITLEN(p)), isfeasible (p, ptr));
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
  uint16_t i, j;
  uint64_t iter;
  
  for (i = 0; i < p->graph->nvert; i++) {
    if (p->graph->vtable[i] == v) {
      if (!(chrom[i / 64] & (1 << (i%64))))
        return 1;
      return 0;
    }
  }
  return 0;
}

float sumweights (pool_s *p, uint64_t *chrom)
{
  uint8_t pos;
  uint16_t i, j, csize;
  uint64_t iter;
  uint64_t *ptr;
  vertex_s *v;
  float weight;
  
  for (weight = 0, i = 0, ptr = chrom; i < p->chromsize; i++, ptr++) {
    csize = (i == p->chromsize-1) ? p->remain : 64;
    for (iter = *ptr, pos = 0; pos <= csize; iter &= ~(1 << pos)) {
      pos = ffsl(iter);
      if (!pos)
        break;
      --pos;
      v = p->graph->vtable[i*64 + pos];
      for (j = 0; j < v->nedges; j++) {
        if (iscut(p,ptr, (v->edges[j]->v1 == v) ? v->edges[j]->v2 : v->edges[j]->v1))
          weight += v->edges[j]->weight;
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
  for (i = 0; i < _POOLSIZE; i++) {
    printf("Weight: %f, %d, %d\n",sumweights (p, ptr), n, isfeasible(p, ptr));
    ptr += n;
  }
}

float getfitness (pool_s *p, uint64_t *chrom)
{
  int setcount, differ;
  
  setcount = countdigits (p, chrom);
  differ = abs(2*setcount-_GET_CHBITLEN(p));
  return sumweights (p, chrom) + (differ<<8);
}

int prcmp (roulette_s *a, roulette_s *b)
{
  if (a->prob < b->prob)
    return -1;
  if (a->prob > b->prob)
    return 1;
  return 0;
}

float computeprob (pool_s *p)
{
  uint16_t i, n;
  uint64_t *ptr;
  float sum;
  roulette_s roul[_POOLSIZE];
  
  for (sum = 0, n = p->chromsize, ptr = p->popul, i = 0; i < _POOLSIZE; i++, ptr += n) {
    roul[i].prob = getfitness (p, ptr);
    roul[i].ptr = ptr;
    sum += roul[i].prob;
  }
  for (i = 0; i < _POOLSIZE; i++) {
    p->rbuf[i].prob = sum / roul[i].prob;
    p->rbuf[i].ptr = roul[i].ptr;
  }
  qsort (p->rbuf, _POOLSIZE, sizeof(roulette_s), (int (*)(const void *, const void *))prcmp);
  return sum;
}

int isfeasible (pool_s *p, uint64_t *chrom)
{
  return (2*countdigits (p, chrom) == _GET_CHBITLEN(p));
}

void singlepoint_cr (pool_s *p, uint64_t *p1, uint64_t *p2)
{
  uint16_t  point,
            remain,
            index, i,
            bitlen;
  uint64_t  tmp1, tmp2, mask;
  
  if ((p1 == p->rbuf[49].ptr && isfeasible(p,p1)) || (p2 == p->rbuf[49].ptr && isfeasible(p, p2)))
    return;
  bitlen = _GET_CHBITLEN(p);
  point = rand() % bitlen;
  index = point % 64;
  tmp1 = *p1;
  tmp2 = *p2;
  for (i = 0, mask = 0; i < index+1; i++)
    mask |= (1 << i);
  *p1 = p2[point / 64] >> (index+1);
  p2[point / 64] = ((tmp1 & mask) << (bitlen-(index+1)));
  *p2 |= tmp1 >> (index+1);
  p1[point / 64] |= ((tmp2 & mask) << (bitlen-(index+1)));
}

void doublepoint_cr (pool_s *p, uint64_t *p1, uint64_t *p2)
{
  
}

void mutate1 (pool_s *p, uint64_t *victim)
{
  
}

#define CBUF_SIZE 32

int run_ge (wgraph_s *g)
{
  pid_t pid;
  pool_s *p;
  uint16_t i1, i2, n;
  int index;
  unsigned char cbuf[CBUF_SIZE];
  gtoken_s *tokens;
  
  if (signal(SIGINT, sigdummy) == SIG_ERR)
    return -1;
  memset (cbuf, 0, sizeof(cbuf));
  pid = fork();
  if (pid) {
    cbuf[CBUF_SIZE-1] = _UEOF;
    p = malloc(sizeof(*p));
    if (!p)
      return -1;
    while (1) {
      index = 0;
      while ((cbuf[index] = (char)getchar()) != '\n') {
        if (index < CBUF_SIZE-2)
          ++index;
      }
      cbuf[index] = '\0';
      cbuf[index+1] = _UEOF;
      tokens = lex_ (cbuf);
      if (!strcmp(cbuf, "new")) {
        printf("Migrating: %d\n", getpid());
      }
      else if (!strcmp(cbuf, "status")) {
        kill (pid, SIGINT);
      }
      else if (!strcmp(cbuf, "show")) {
        
      }
      else if (!(
                  strcmp(cbuf, "exit")  &&
                  strcmp(cbuf, "Q")     &&
                  strcmp(cbuf, "q")     &&
                  strcmp(cbuf, "quit")
                )
      )
      {
        printf("Final:\n");
        kill(pid, SIGINT);
        pause();
        kill(pid, SIGQUIT);
        exit(EXIT_SUCCESS);
      }
      freetokens (tokens);
    }
  }
  else {
    if (signal(SIGINT, siginthandle) == SIG_ERR)
      return -1;
    if (signal(SIGQUIT, sigquithandle) == SIG_ERR)
      return -1;
    p = pool_init (g);
    n = p->chromsize;
    while (1) {
      computeprob (p);
      i1 = rand()%_POOLSIZE;
      while ((i2 = rand()%_POOLSIZE) == i1);
      p->cross (p, p->rbuf[i1*n].ptr, p->rbuf[i2*n].ptr);
      p->mutate (p, p->popul);
      ++p->gen;
    }
  }
  return 0;
}

void siginthandle (int signal)
{
  printf("Status at Generation: %llu\n", pool_->gen);
  printpool(pool_);
  printf("\nMost Fit ( weight = %f ):\n",sumweights(pool_,pool_->rbuf[49].ptr));
  printchrom (pool_, pool_->rbuf[49].ptr);
  if (isfeasible(pool_, pool_->rbuf[49].ptr))
      printf("\nIs Feasible\n");
  else
    printf("\nNot Feasible\n");
  kill(getppid(), SIGINT);
}

void sigquithandle (int signal)
{
  exit(EXIT_SUCCESS);
}

void sigdummy (int signal){}
