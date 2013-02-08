#include "bis.h"
#include <stdlib.h>

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
