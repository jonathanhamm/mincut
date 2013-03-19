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

/*
 Counts the set bits in a long word:
 32-bit count code obtained from:
 http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
 */
#define countbitsLW(lw) ((lw & 0xfff) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f \
+ (((lw & 0xfff000) >> 12) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f \
+ ((lw >> 24) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f

static uint8_t bmap8[256] = {
  0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
};

static pool_s *pool_s_ (uint16_t csize);
static void printlword (uint64_t lword, uint8_t mask);
static float sumweights (pool_s *p, uint64_t *chrom);
static float getfitness (pool_s *p, uint64_t *chrom);
static float computeprob (pool_s *p);
static int isfeasible (pool_s *p, uint64_t *chrom);

pool_s *pool_s_ (uint16_t csize)
{
  int i;
  pool_s *pool;
  
  pool = calloc(1, sizeof(*pool) + _POOLSIZE*_CQWORDSIZE(csize)*8);
  if (!pool)
    goto err_;
  pool->chromsize = (csize >> 6) + (csize % 64 != 0);
  pool->remain = csize % 64;
  for (i = 0; i < pool->remain; i++)
    pool->cmask |= (1 << i);
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

void printpool (pool_s *p)
{
  uint16_t i, n, index;
  uint64_t *ptr;
  
  n = p->chromsize;
  ptr = p->popul;
  for (index = 0; index < _POOLSIZE; index++) {
    printf("%d:", index);
    for (i = 0; i < n; i++) {
      if (i == n-1)
        printlword (ptr[i], p->remain);
      else
        printlword (ptr[i], 64);
    }
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
  return sumweights (p, chrom) + (differ<<4);
}

float computeprob (pool_s *p)
{
  uint16_t i, n;
  uint64_t *ptr;
  float sum;
  float fitbuf[_POOLSIZE];
  
  for (sum = 0, n = p->chromsize, ptr = p->popul, i = 0; i < _POOLSIZE; i++, ptr += n) {
    fitbuf[i] = getfitness (p, ptr);
    sum += fitbuf[i];
  }
  for (i = 0; i < _POOLSIZE; i++)
    p->probuf[i] = sum / fitbuf[i];
  return sum;
}

int isfeasible (pool_s *p, uint64_t *chrom)
{
  return (2*countdigits (p, chrom) == _GET_CHBITLEN(p));
}

