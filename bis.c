#include "bis.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define __UEOF (unsigned char)EOF
#define __INITBUFSIZE 256
#define __MAXLEXLEN 16

/*TTYPES*/
#define _ID 1
#define _OPENBRACE 2
#define _CLOSEBRACE 3
#define _COMMA 4

typedef struct gtoken_s gtoken_s;

struct gtoken_s
{
    unsigned short type;
    unsigned char lexeme[__MAXLEXLEN + 1];
    gtoken_s *prev;
    gtoken_s *next; 
};

static gtoken_s *stream_;

static gtoken_s *gtoken_s_ (gtoken_s *node, unsigned char *lexeme, unsigned short type);
static gtoken_s *lex_ (unsigned char *buf);

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


/*  Graph Parsing Routines */

/*  Graph Grammar:
 *  <graph> => <nodelist> <edgelist> EOF
 *  <nodelist> => N={n <nodeparam>}
 *  <nodeparam> => ,n <nodeparam> | E
 *  <edgelist> => E={<e> <edgeparam> }
 *  <edgeparam> => ,<e> <edgeparam> | E
 *  <e> => {n,n,real}
 *
 *  Regex for node/edge definitions:
 *  <N> : (a...Z)+(a...z+0...9)*
 *  <e> : (0...9)+<optional_fract>
 *  <optional_fract> : (0...9)?
 */
gnode_s *gparse (unsigned char *buf)
{
    
    if (!lex_ (buf))
        return NULL;
}

gtoken_s *gtoken_s_ (gtoken_s *node, unsigned char *lexeme, unsigned short type)
{
    gtoken_s *tok;
    
    tok = calloc(1,sizeof(*tok));
    if (!tok)
        return NULL;
    tok->type = type;
    strcpy(tok->lexeme, lexeme);
    if (!node)
        stream_ = tok;
     else {
        tok->prev = node;
        node->next = tok;
    }
    return tok;
}

gtoken_s *lex_ (unsigned char *buf)
{
    unsigned char   backup;
    unsigned char   *bckptr;
    gtoken_s        *curr;
    
    for (curr = NULL, bckptr = buf; *buf != __UEOF;) {
        switch (*buf) {
            case ' ':
            case '\n':
            case '\r':
            case '\t':
            case '\0':
            case '\v':
                continue;
            case ',':
                curr = gtoken_s_ (curr, ",", _COMMA);
                break;
            default:
                if ((*buf > 'A' && *buf < 'Z') || (*buf > 'a' && *buf < 'z')) {
                    for (bckptr = buf++; (*buf >= 'A' && *buf <= 'Z') || (*buf >= 'a' && *buf <= 'z')
                         || (*buf >= '0' && *buf <='9'); buf++) {
                        if (buf - bckptr == __MAXLEXLEN) {
                            printf("Too Long Lexeme: %16s", bckptr);
                            goto err_;
                        }
                    }
                    backup = *buf;
                    *buf = '\0';
                    curr = gtoken_s_ (curr, bckptr, _ID);
                    *buf = backup;
                    bckptr = buf;
                } else if (*buf > '0' && *buf < '9') {
                    for (buf++; (*buf > '0' && *buf < '9'); buf++) {
                        if (buf - bckptr == __MAXLEXLEN) {
                            printf("Too Long Lexeme: %16s", bckptr);
                            goto err_;
                        }
                    }
                    if (*buf == '.') {
                        for (buf++; (*buf > '0' && *buf < '9'); buf++) {
                            if (buf - bckptr == __MAXLEXLEN) {
                                printf("Too Long Lexeme: %16s", bckptr);
                                goto err_;
                            }
                        }
                    }
                    backup = *buf;
                    *buf = '\0';
                    curr = gtoken_s_ (curr, bckptr, _ID);
                    *buf = backup;
                    bckptr = buf;
                }
        }
    }
err_:
    return NULL;
}

