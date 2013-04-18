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
static void sigmutatehandle (int signal);

static pool_s *pool_s_ (uint16_t csize);
static void printlword (uint64_t lword, uint8_t mask);
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

pool_s *pool_s_ (uint16_t csize)
{
  int i;
  pool_s *p;
  
  p = calloc(1, sizeof(*p) + _POOLSIZE*_CQWORDSIZE(csize)*8);
  if (!p)
    goto err_;
  pool_ = p;
  p->chromsize = (csize / 64) + (csize % 64 != 0);
  p->remain = csize % 64;
  if (!(csize % 64) && csize)
    p->cmask = 0xffffffffffffffffllu;
  else {
    for (i = 0; i < p->remain; i++)
        p->cmask |= (1llu << i);
  }
  p->select = roulette_sf;
  p->cross = uniform_cr;
  p->mutate = mutate1;
  printf("QWORD size: %d\n", p->chromsize);
  printf("Remainder: %d\nMask\n", p->remain);
  printlword(p->cmask, 64);
  printf ("\n");
  return p;
  
err_:
  return NULL;
}

pool_s *pool_init (wgraph_s *g)
{
  uint16_t  i, j;
  pool_s    *p;
  uint64_t  *ptr;
  float sum;
  
  p = pool_s_(g->nvert);
  if (!p)
    return NULL;
  p->gen = 0;
  ptr = p->popul;
  for (i = 0; i < _POOLSIZE; i++) {
    for (j = 0; j < p->chromsize; j++, ptr++) {
      ((uint16_t *)ptr)[0] = (uint16_t)rand();
      ((uint16_t *)ptr)[1] = (uint16_t)rand();
      ((uint16_t *)ptr)[2] = (uint16_t)rand();
      ((uint16_t *)ptr)[3] = (uint16_t)rand();
    }
     *(ptr-1) &= p->cmask;
  }
  p->graph = g;
  p->mutateprob = _INITMUTATIONPROB;
  p->bitlen = ((p->remain) ? ((p->chromsize * 64) - (64 - p->remain)) : p->chromsize*64);
  for (sum = 0, ptr = p->popul, i = 0; i < _POOLSIZE; i++, ptr += p->chromsize) {
    p->rbuf[i].fitness = getfitness (p, ptr);
    p->rbuf[i].ptr = ptr;
    sum += p->rbuf[i].fitness;
  }
  for (i = 0, p->accum = 0; i < _POOLSIZE; i++) {
    p->rbuf[i].prob = sum / p->rbuf[i].fitness;
  }
  qsort (p->rbuf, _POOLSIZE, sizeof(roulette_s), (int (*)(const void *, const void *))prcmp);
  for (i = 0; i < _POOLSIZE; i++) {
    p->accum += p->rbuf[i].prob;
    p->rbuf[i].cummulative = (int)p->rbuf[i].prob;
    if (i)
      p->rbuf[i].cummulative += p->rbuf[i-1].cummulative;
  }
  p->fitsum = sum;
  return p;
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

int bsearch_r (roulette_s *roul, uint32_t key)
{
	int mid, low, high;

  low = 0;
  high = _POOLSIZE-1;
  for (mid = low+(high-low)/2; high >= low; mid = low+(high-low)/2) {
    if (key < roul[mid].cummulative)
      high = mid - 1;
    else if (key > roul[mid].cummulative)
      low = mid + 1;
    else
      return mid;
  }
  return mid;
}

void roulette_sf (pool_s *p, selected_s *parents)
{
  int i;
  uint32_t i1, i2;
  
  for (i = 0; i < _NSELECT; i++) {
    i1 = bsearch_r (p->rbuf, rand() % p->accum);
    i2 = bsearch_r (p->rbuf, rand() % p->accum);
    if (i2 == i1)
      i2 = (i2 - 1) % _POOLSIZE;
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
  uint64_t iter;
  uint64_t *ptr;
  vertex_s *v;
  edge_s **edges;
  float weight;

  for (weight = 0, i = 0, ptr = chrom; i < p->chromsize; i++, ptr++) {
    if (i == p->chromsize-1 && p->remain)
      csize = p->remain;
    else
      csize = 63;
    for (iter = *ptr, pos = 0; pos < csize; iter &= ~(1llu << pos)) {
      pos = ffsl(iter);
      if (!pos)
        break;
      --pos;
      v = p->graph->vtable[/* i*64 */ (i << 6) + pos];
      nedges = v->nedges;
      edges = v->edges;
      for (j = 0; j < nedges; j++) {
        if (iscut(p, ptr, (edges[j]->v1 == v) ? edges[j]->v2 : edges[j]->v1))
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
  return sumweights (p, chrom) + (differ << 4);
}

void computeprob (pool_s *p)
{
  uint16_t i, n;
  uint64_t *ptr;
  float sum, fc1, fc2;
  
  sum = p->fitsum;
  for (i = 0, p->accum = 0; i < _POOLSIZE; i++)
    p->rbuf[i].prob = sum / p->rbuf[i].fitness;
  qsort (p->rbuf, _POOLSIZE, sizeof(roulette_s), (int (*)(const void *, const void *))prcmp);
  for (i = 0; i < _POOLSIZE; i++) {
    p->rbuf[i].cummulative = (int)p->rbuf[i].prob;
    if (i)
      p->rbuf[i].cummulative += p->rbuf[i-1].cummulative;
    p->accum += p->rbuf[i].prob;
  }
  if (isfeasible(p, p->rbuf[_POOLSIZE-1].ptr))
    p->bestfeasible = p->rbuf[_POOLSIZE-1].ptr;
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
  
  if ((p1 == p->rbuf[_POOLSIZE-1].ptr && isfeasible(p,p1)) || (p2 == p->rbuf[_POOLSIZE-1].ptr && isfeasible(p, p2)))
    return;
  bitlen = _GET_CHBITLEN(p);
  point = rand() % bitlen;
  index = point % 64;
  tmp1 = *p1;
  tmp2 = *p2;
  for (i = 0, mask = 0; i < index+1; i++)
    mask |= (1llu << i);
  *p1 = p2[point / 64] >> (index+1);
  p2[point / 64] = ((tmp1 & mask) << (bitlen-(index+1)));
  *p2 |= tmp1 >> (index+1);
  p1[point / 64] |= ((tmp2 & mask) << (bitlen-(index+1)));
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
    rp2 = &p->rbuf[1%_POOLSIZE];
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
  
  n = _GET_CHBITLEN(p) >> 4;
  if (!isfeasible(p, victim))
    ++n;
  for (i = 0; i < n; i++) {
    index = rand() % _GET_CHBITLEN(p);
    victim[index / 64] ^= (1llu << (index % 64));
  }
}

void mutate2 (pool_s *p, uint64_t *victim)
{
  uint16_t index, i;
  
  for (i = 0; i < 10; i++) {
    index = rand() % _GET_CHBITLEN(p);
    victim[index / 64] ^= (1llu << (index ));
  }
}

#define CBUF_SIZE 32


void pSIGALRM (int signal) {}

int run_ge (wgraph_s *g)
{
  pool_s *p;
  int i;
  uint16_t n, index;
  selected_s parents;
  unsigned char cbuf[CBUF_SIZE];
  
  if (signal(SIGINT, sigdummy) == SIG_ERR)
    return -1;
  if (signal(SIGALRM, pSIGALRM) == SIG_ERR)
    return -1;
  memset (cbuf, 0, sizeof(cbuf));
  p = pool_init (g);
  if (!p)
    return -1;
  if (pipe (pipe_) == -1)
    return -1;
  pid_ = fork();
  if (pid_) {
    cbuf[CBUF_SIZE-1] = _UEOF;
    while (1) {
      index = 0;
      while ((cbuf[index] = (char)getchar()) != '\n') {
        if (index < CBUF_SIZE-1)
          ++index;
        else
          printf ("Command Length %d Exceeded\n", CBUF_SIZE-1);
      }
      cbuf[index] = _UEOF;
      write(pipe_[1], cbuf, CBUF_SIZE);
      kill (pid_, SIGUSR1);
      pause ();
    }
  }
  else {
    if (signal(SIGUSR1, siginthandle) == SIG_ERR)
      return -1;
    if (signal(SIGINT, pSIGALRM) == SIG_ERR)
      return -1;
    n = p->chromsize;
    p->start = clock();
    while (1) {
      p->select (p, &parents);
      for (i = 0; i < _NSELECT; i++)
        p->cross (p, parents.couples[i].p1, parents.couples[i].p2);
      computeprob(p);
      ++p->gen;
    }
  }
  return 0;
}

void printstatus (void)
{
  int i;
  
  printf("Status at Generation: %llu\n", pool_->gen);
  printpool(pool_);
  printf("\nMost Fit ( weight = %f ):\n", pool_->rbuf[_POOLSIZE-1].fitness);
  printchrom (pool_, pool_->rbuf[_POOLSIZE-1].ptr);
  if (isfeasible(pool_, pool_->rbuf[_POOLSIZE-1].ptr))
    printf("\nIs Feasible\n");
  else {
    printf("\nNot Feasible\n");
    if (pool_->bestfeasible) {
      printf("\nMost Fit Feasible ( weight = %f ):\n", getfitness(pool_, pool_->bestfeasible));
          printchrom (pool_, pool_->rbuf[i].ptr);
    }
    else
      printf("No Feasibles Found\n");
  }
  printf("Time Elapse: %llu\n", (uint64_t)(clock() - pool_->start));

}

void siginthandle (int signal)
{
  int tmpint;
  unsigned char cbuf[CBUF_SIZE];
  gtoken_s *head, *tokens;
  
  printf("called\n");
  read(pipe_[0], cbuf, CBUF_SIZE);
  head = lex_ (cbuf);
  tokens = head;
  if (!tokens) {
    printf ("Illegal Symbols Used\n");
    goto exit_;
  }
  else if (!strcmp(tokens->lexeme, "status")) {
    printstatus ();
  }
  else if (!strcmp(tokens->lexeme, "set")) {
    tokens = tokens->next;
    if (!strcmp(tokens->lexeme, "mutate")) {
      tokens = tokens->next;
      if (!strcmp(tokens->lexeme, "prob")) {
        tokens = tokens->next;
        if (tokens->type == _NUM) {
          tmpint = atoi (tokens->lexeme);
          if (tmpint > 100)
            tmpint = 100;
          else if (tmpint < 0)
            tmpint = 0;
          printf("Mutation Probability Set at: %d%%\n", tmpint);
          pool_->mutateprob = tmpint;
        }
      }
      else if (!strcmp(tokens->lexeme, "op")) {
        tokens = tokens->next;
        if (tokens->type == _NUM) {
          tmpint = atoi (tokens->lexeme);
          if (tmpint == 1) {
            printf ("Set to Mutate Function 1\n");
            pool_->mutate = mutate1;
          }
          else if (tmpint == 2) {
            printf ("Set to Mutate Function 2\n");
            pool_->mutate = mutate2;
          }
          else
            printf ("Expected 1 or 2, but got %s\n", tokens->lexeme);
        }
        else 
          printf ("Expected Number, but got %s\n", tokens->lexeme);
      }
      else
        printf ("Expected 'op' or 'prob' but got %s\n", tokens->lexeme);
    }
  }
  else if (!strcmp(tokens->lexeme, "get")) {
    tokens = tokens->next;
    if (!strcmp(tokens->lexeme, "mutate")) {
      tokens = tokens->next;
      if (!strcmp(tokens->lexeme, "prob"))
        printf ("Mutation Probability: %d%%\n", pool_->mutateprob);
      else if (!strcmp(tokens->lexeme, "op")) {
        if (pool_->mutate == mutate1)
          printf("Mutation Operator 1\n");
        else
          printf("Mutation Operator 2\n");
      }
      else
        printf ("Expected 'probability' or 'operator', but got %s\n", tokens->lexeme);
    }
    else
      printf ("Expected Attribute Name, but got: %s\n", tokens->lexeme);
  }
  else if (!(
             strcmp(tokens->lexeme, "exit")  &&
             strcmp(tokens->lexeme, "Q")     &&
             strcmp(tokens->lexeme, "q")     &&
             strcmp(tokens->lexeme, "quit")
             )
           )
  {
    printf("Final:\n");
    printstatus ();
    kill(getppid(), SIGQUIT);
    exit(EXIT_SUCCESS);
  }
  else {
    printf("'%s' not recognized.\n",tokens->lexeme);
  }
  
exit_:
  freetokens (head);
  kill(getppid(), SIGALRM);
}

void sigdummy (int signal)
{
  unsigned char cbuf[CBUF_SIZE];
  
  strcpy(cbuf, "status");
  cbuf[6] = _UEOF;
  write(pipe_[1], cbuf, CBUF_SIZE);
  kill (pid_, SIGUSR1);
  pause();
}

void printsolution (pool_s *p, int index)
{
  
}

inline void printfittest (pool_s *p)
{
  printsolution (p, _POOLSIZE-1);
}
