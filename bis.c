/*
* Main GA Engine
* Author: Jonathan Hamm
*/

#include "gparse.h"
#include "bis.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

struct chromosome_s
{
  uint64_t mask;
  uint64_t data;
};

pool_s *pool_s_ (uint16_t csize);
void printlword (uint64_t lword, uint8_t mask);

pool_s *pool_s_ (uint16_t csize)
{
  pool_s *pool;
  
  pool = calloc(1, sizeof(*pool) + _POOLSIZE*_CQWORDSIZE(csize)*8);
  if (!pool)
    goto err_;
  pool->chromsize = (csize >> 6) + (csize % 64 != 0);
  pool->cmask = csize % 64;
  return pool;
  
err_:
  return NULL;
}

pool_s *pool_init (wgraph_s *g)
{
  int i, n;
  pool_s *pool;
  
  printf("nedges: %d\n", g->nedges);
  pool = pool_s_(g->nedges);
  if (!pool)
    return NULL;
  n = pool->chromsize*4*_POOLSIZE;
  for (i = 0; i < n; i++)
    ((uint16_t *)&pool->popul)[i] = (uint16_t)rand();  /*rand() returns a 4 byte integer, but it ignores the signed bit*/
  return pool;
}

void printlword (uint64_t lword, uint8_t end)
{
  uint8_t i;
  
  if (!end)
    end = 64;
  for (i = 0; i < end; i++)
    printf("%d",(int)(((lword << i)) >> 63));
}

void printpool (pool_s *p)
{
  uint16_t i, n, index;
  uint64_t *ptr;
  
  n = p->chromsize;
  for (index = 0; index < _POOLSIZE; index++) {
    ptr = p->popul + index * n;
    for (i = 0; i < n; i++) {
      if (i == n-1)
        printlword (ptr[i], p->cmask);
      else
        printlword (ptr[i], 64);
    }
    printf("\n");
  }
}

