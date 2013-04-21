#include <stdio.h>

#define NCHROM 512

int main(int argc, const char * argv[])
{
  int i, j;

  printf ("V = \n{");
  for (i = 0; i < NCHROM-1; i++) {
    if (!(i%16))
      printf("\n\t");
    printf("v%03d, ",i);
  }
  printf("v%03d\n}\n\nE = {\n", NCHROM-1);
  for (i = 0; i < NCHROM/2; i++)
    printf("\t{v%03d, v%03d, 1},\n", i, i+NCHROM/2);
  for (i = 0; i < NCHROM/2; i++) {
    for (j = i+1; j < NCHROM/2; j++)
      printf("\t{v%03d, v%03d, 1.01},\n", i, j);
  }
  for (; i < NCHROM-2; i++) {
    for (j = i+1; j < NCHROM; j++)
      printf("\t{v%03d, v%03d, 1.01},\n", i, j);
  }
  printf("\t{v%02d, v%03d, 1.01}\n}", i, j-1);
  return 0;
}

