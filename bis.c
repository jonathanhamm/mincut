#include "bis.h"
#include <stdlib.h>
#include <stdio.h>

#define __UEOF (unsigned char)EOF
#define __INITBUFSIZE 256

pool_s *pool_s_ (uint16_t csize) 
{
  pool_s *pool;

  pool = malloc(sizeof(*pool));
  if (!pool)
    goto err_;

err_:
  if (pool)
    free(pool);
  return NULL;
}

unsigned char *read_gfile (const char *fname) 
{
  FILE *f;
  size_t offset,
         bsize;
  unsigned char *buf;
              

  f = fopen(fname,"r");
  if(!f) {
    printf("Error Opening: %s\n",fname);
    return NULL;
  }
  buf = malloc(__INITBUFSIZE);
  if (!buf)
    goto err_;
  for (bsize = __INITBUFSIZE, offset = 0; (buf[offset] = (unsigned char)fgetc(f)) != __UEOF; offset++) {
    if (offset == bsize) {
      bsize += __INITBUFSIZE;
      buf = realloc (buf, bsize);
      if (!buf)
        goto err_;
    } 
  }
  /*truncate buffer to EOF*/
  if (offset < bsize)
    buf = realloc (buf, offset+1);
  return buf;

err_:
  printf("Heap Allocation Error\n");
  fclose(f);
  return NULL;
}
