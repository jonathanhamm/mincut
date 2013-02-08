#include "bis.h"
#include <stdio.h>

int main (int argc, char **argv) 
{
  unsigned char *tmp;

  if (argc != 2) {
    printf("Usage: %s <filename>", argv[0]);
    return 1;
  }
  tmp = read_gfile (argv[1]);
  if (!tmp)
    return 2;
  return 0;
}
