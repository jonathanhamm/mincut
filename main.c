#include "gparse.h"
#include "bis.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

int main (int argc, char **argv) 
{
  clock_t t;
  
  if (argc != 2) {
    printf("Usage: ./ge <filename>");
    exit (EXIT_FAILURE);
  }
  //run_ge (gparse(argv[1]));
  run_simanneal (gparse(argv[1]));
  exit (EXIT_SUCCESS);
}
