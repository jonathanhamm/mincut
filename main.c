#include "gparse.h"
#include "bis.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

int main (int argc, char **argv) 
{
  pool_s *p;
  clock_t t;
  int i;
  
  if (argc != 2) {
    printf("Usage: %s <filename>", argv[0]);
    return 1;
  }
  p = pool_init (gparse (argv[1]));
  printpool (p);
  int x = 0;
  clock_t fin;
  clock_t sum = 0;
  for (i = 0; i < 1000000; i++) {
    t = clock();
    //x += countdigits(p,rand()%50);
    fin = clock() - t;
    sum+=fin;
  }
  printf("time: %lu, %f, %d\n", sum, sum/1000000.0, x);
  printweights(p);
  return 0;
}
