#include "parse.h"
#include "bis.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

int main (int argc, char **argv) 
{  
  if (argc != 2) {
    if (argc == 3) {
      if (!strcmp(argv[2], "ge"))
        run_ge (gparse(argv[1]));
      else if (!strcmp(argv[2], "sa"))
        run_simanneal (gparse(argv[1]), SIMULATED_ANNEALING);
      else if (!strcmp(argv[2], "hc"))
        run_simanneal (gparse(argv[1]), HILL_CLIMBNG);
      else {
        printf("Usage: ./ge <filename> <algorithm>");
        exit (EXIT_FAILURE);
      }
    }
    else {
      printf("Usage: ./ge <filename> <algorithm>");
      exit (EXIT_FAILURE);
    }
  }
  else
    run_ge (gparse(argv[1]));
  exit (EXIT_SUCCESS);
}
