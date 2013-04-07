#include "gparse.h"
#include "bis.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

int main (int argc, char **argv) 
{
  clock_t t;
  
  if (argc != 2) {
    printf("Usage: %s <filename>", argv[0]);
    return 1;
  }
  run_ge (gparse(argv[1]));
  return 0;
}
