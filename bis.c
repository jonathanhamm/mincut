#include "bis.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define _UEOF (unsigned char)EOF
#define _INITBUFSIZE 256
#define _MAXLEXLEN 16

/*TTYPES*/
#define _ID         1
#define _OPENBRACE  2
#define _CLOSEBRACE 3
#define _COMMA      4
#define _EQU        5
#define _EOF        6

#define __GTNEXT(str) (str = str->next)
#define __GTPREV(str) (str = str->prev)

typedef struct gtoken_s gtoken_s;

struct gtoken_s
{
    unsigned short type;
    unsigned char lexeme[_MAXLEXLEN + 1];
    gtoken_s *prev;
    gtoken_s *next; 
};

static gtoken_s *stream_;

static gtoken_s *gtoken_s_ (gtoken_s *node, unsigned char *lexeme, unsigned short type);
static gtoken_s *lex_ (unsigned char *buf);
static int parse_ (void);

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
  buf = malloc(_INITBUFSIZE);
  if (!buf)
    goto err_;
  for (bsize = _INITBUFSIZE, offset = 0; (buf[offset] = (unsigned char)fgetc(f)) != _UEOF; offset++) {
    if (offset == bsize) {
      bsize += _INITBUFSIZE;
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

gnode_s *gparse (unsigned char *buf)
{
    
    if (!lex_ (buf))
        return NULL;
}

/* 
*  Regex for node/edge definition lexemes:
*  n: (a...Z)+(a...z+0...9)*
*  real: (0...9)+<optional_fract>
*  <optional_fract>: (dot 0...9)?
*/
gtoken_s *lex_ (unsigned char *buf)
{
    unsigned char   backup;
    unsigned char   *bckptr;
    gtoken_s        *curr;
    
    for (curr = NULL, bckptr = buf; *buf != _UEOF;) {
        switch (*buf) {
            case ',':
                curr = gtoken_s_ (curr, ",", _COMMA);
                buf++;
                break;
            case '=':
                curr = gtoken_s_ (curr, "=", _EQU);
                buf++;
                break;
            case '{':
                curr = gtoken_s_ (curr, "{", _OPENBRACE);
                buf++;
                break;
            case '}':
                buf++;
                curr = gtoken_s_ (curr, "}", _CLOSEBRACE);
                break;
            default:
                if (*buf <= ' ')
                    buf++;
                 else if ((*buf >= 'A' && *buf <= 'Z') || (*buf >= 'a' && *buf <= 'z')) {
                    for (bckptr = buf, ++buf; (*buf >= 'A' && *buf <= 'Z') || (*buf >= 'a' && *buf <= 'z')
                         || (*buf >= '0' && *buf <= '9'); buf++) {
                        if (buf - bckptr == _MAXLEXLEN) {
                            printf("Too Long Lexeme: %16s", bckptr);
                            goto err_;
                        }
                    }
                    backup = *buf;
                    *buf = '\0';
                    curr = gtoken_s_ (curr, bckptr, _ID);
                    *buf = backup;
                    bckptr = buf;
                } else if (*buf >= '0' && *buf <= '9') {
                    for (bckptr = buf, buf++; (*buf >= '0' && *buf <= '9'); buf++) {
                        if (buf - bckptr == _MAXLEXLEN) {
                            printf("Too Long Lexeme: %16s", bckptr);
                            goto err_;
                        }
                    }
                    if (*buf == '.') {
                        for (buf++; (*buf >='0' && *buf <= '9'); buf++) {
                            if (buf - bckptr == _MAXLEXLEN) {
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
                } else {
                    printf("Symbol Error\n");
                    goto err_;
                }
                break;
        }
    }
    curr = gtoken_s_ (curr, "\xFF", _EOF);
    return stream_;
err_:
    return NULL;
}

/*  Graph Grammar:
 *  <graph> => <nodelist> <edgelist> EOF
 *  <nodelist> => N={n <nodeparam>}
 *  <nodeparam> => ,n <nodeparam> | E
 *  <edgelist> => E={<e> <edgeparam> }
 *  <edgeparam> => ,<e> <edgeparam> | E
 *  <e> => {n,n,real}
 */
static int parse_ (void)
{
    __GTNEXT(str)
}