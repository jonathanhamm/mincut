/*
 main.c 
 Author: Jonathan Hamm
 
 Description:
 
 This simply has the main function of the progran. It 
 accepts commandline arguments which specify what file to
 open (containing the graph data), and what algorithm to use. 
 It then calls the appropriate functions. 
 */

#include "parse.h"
#include "bis.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

int main (int argc, char **argv)
{
    if (argc != 2) {
        if (argc == 4) {
            if (!strcmp(argv[2], "-a")) {
                if (!strcmp(argv[3], "ge"))
                    run_ge (gparse(argv[1]));
                else if (!strcmp(argv[3], "sa"))
                    run_simanneal (gparse(argv[1]), SIMULATED_ANNEALING);
                else if (!strcmp(argv[3], "hc"))
                    run_simanneal (gparse(argv[1]), HILL_CLIMBNG);
                else {
                    printf("Usage: ./ge <filename> -a <algorithm>");
                    exit (EXIT_FAILURE);
                }
            }
        }
        else {
            printf("Usage: ./ge <filename> -a <algorithm>");
            exit (EXIT_FAILURE);
        }
    }
    else
        run_ge (gparse(argv[1]));
    exit (EXIT_SUCCESS);
}
