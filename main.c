#include "gparse.h"
#include "bis.h"
#include <stdio.h>

int main (int argc, char **argv) 
{
  pool_s *p;

  if (argc != 2) {
    printf("Usage: %s <filename>", argv[0]);
    return 1;
  }
  p = pool_init (gparse (argv[1]));
  printpool (p);
  return 0;
}
