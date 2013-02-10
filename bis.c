/*
 * Main GA Engine
 * Author: Jonathan Hamm
 */

#include "bis.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define _VTABLE_SIZE 19
#define _CHUNK_SIZE 8
#define _UEOF (unsigned char)EOF
#define _INITBUFSIZE 256
#define _MAXLEXLEN 15

/*TTYPES*/
#define _ID         1
#define _NUM        2
#define _OPENBRACE  3
#define _CLOSEBRACE 4
#define _COMMA      5
#define _EQU        6
#define _EOF        7

#define __GTNEXT() (stream_ = stream_->next)
#define __GTPREV() (stream_ = stream_->prev)

typedef struct gtoken_s gtoken_s;
typedef struct wgraph_s wgraph_s;
typedef struct vrecord_s vrecord_s;
typedef struct vertex_s vertex_s;
typedef struct vchain_s vchain_s;
typedef struct edge_s edge_s;

struct gtoken_s
{
    unsigned short type;
    unsigned char lexeme[_MAXLEXLEN + 1];
    gtoken_s *prev;
    gtoken_s *next; 
};

struct vrecord_s
{
    vertex_s *v;
    union {
        unsigned long isoccupied;
        vchain_s *chain;
    };
};

struct vchain_s
{
    uint8_t mem;
    vrecord_s chunk[_CHUNK_SIZE];
    vchain_s *next;
};

struct vertex_s
{
    unsigned char name[_MAXLEXLEN + 1];
    short nedges;
    edge_s *edges;
};

struct edge_s
{
    double weight;
    vertex_s *v1;
    vertex_s *v2;
};

struct wgraph_s
{
    short nvert;
    vrecord_s vtable[_VTABLE_SIZE];
};


static gtoken_s *stream_;

/*tokenizing routines for graph data*/
static gtoken_s *gtoken_s_ (gtoken_s *node, unsigned char *lexeme, unsigned short type);
static gtoken_s *lex_ (unsigned char *buf);

/*graph parsing routines*/
static int parse_ (void);
static void pgraph_ (void);
static void pnodelist_ (void);
static void pnodeparam_ (void);
static void pedgelist_ (void);
static void pedgeparam_ (void);
static void e_ (void);

/*graph data structure routines*/
static inline wgraph_s *wgraph_s_ (void);
static int insert_vertex (wgraph_s *graph, vertex_s *v);
static vchain_s *vchain_s_ (void);
static int chain_insert (vchain_s *chunk, vertex_s *v);

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
    parse_();
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
                            printf("Too Long ID: %15s", bckptr);
                            goto err_;
                        }
                    }
                    backup = *buf;
                    *buf = '\0';
                    curr = gtoken_s_ (curr, bckptr, _ID);
                    *buf = backup;
                } else if (*buf >= '0' && *buf <= '9') {
                    for (bckptr = buf, buf++; (*buf >= '0' && *buf <= '9'); buf++) {
                        if (buf - bckptr == _MAXLEXLEN) {
                            printf("Too Long ID: %15s", bckptr);
                            goto err_;
                        }
                    }
                    if (*buf == '.') {
                        for (buf++; (*buf >='0' && *buf <= '9'); buf++) {
                            if (buf - bckptr == _MAXLEXLEN) {
                                printf("Too Long ID: %15s", bckptr);
                                goto err_;
                            }
                        }
                    }
                    backup = *buf;
                    *buf = '\0';
                    curr = gtoken_s_ (curr, bckptr, _NUM);
                    *buf = backup;
                } else {
                    printf("Symbol Error\n");
                    goto err_;
                }
                break;
        }
    }
    curr = gtoken_s_ (curr, "\xff", _EOF);
    return stream_;
err_:
    return NULL;
}

/*  Graph Grammar:
 *  <graph> => <nodelist> <edgelist> EOF
 *  <nodelist> => V={v <nodeparam>}
 *  <nodeparam> => ,v <nodeparam> | E
 *  <edgelist> => E={<e> <edgeparam> }
 *  <edgeparam> => ,<e> <edgeparam> | E
 *  <e> => {n,n,real}
 */
static int parse_ (void)
{
    pgraph_();
}

void pgraph_ (void)
{
    pnodelist_();
    pedgelist_();
    if (__GTNEXT()->type == _EOF)
        printf("Parse Success!\n");
}

void pnodelist_ (void)
{
    if (!strcmp(stream_->lexeme, "V"))
    if (__GTNEXT()->type == _EQU)
    if (__GTNEXT()->type == _OPENBRACE)
    if (__GTNEXT()->type == _ID)
        pnodeparam_();
    if (stream_->type == _CLOSEBRACE)
        return;
}

void pnodeparam_ (void)
{
    if (__GTNEXT()->type == _COMMA)
    if (__GTNEXT()->type == _ID)
        pnodeparam_();
}

void pedgelist_ (void)
{
    if (!strcmp(__GTNEXT()->lexeme, "E"))
    if (__GTNEXT()->type == _EQU)
    if (__GTNEXT()->type == _OPENBRACE)
        e_();
    pedgeparam_();
    if (stream_->type == _CLOSEBRACE)
        return;
}

void pedgeparam_ (void)
{
    if (__GTNEXT()->type == _COMMA) {
        e_();
        pedgeparam_();
    }
}

void e_ (void)
{
    if (__GTNEXT()->type == _OPENBRACE)
    if (__GTNEXT()->type == _ID)
    if (__GTNEXT()->type == _COMMA)
    if (__GTNEXT()->type == _ID)
    if (__GTNEXT()->type == _COMMA)
    if (__GTNEXT()->type == _NUM)
    if (__GTNEXT()->type == _CLOSEBRACE)
        return;
}

/*graph data structure routines*/

inline wgraph_s *wgraph_s_ (void)
{
    return calloc(1,sizeof(wgraph_s));
}

int insert_vertex (wgraph_s *graph, vertex_s *v)
{
    uint8_t index;
    vrecord_s *rec;
    
    rec = &graph->vtable[*(uint64_t *)v->name % _VTABLE_SIZE];
    if (!rec->isoccupied) {
        rec->v = v;
        rec->isoccupied = 1;
    } else if (rec->isoccupied == 1) {
        rec->chain = calloc(1,sizeof(*rec->chain));
        
      //  static int chain_insert (vchain_s *chunk, vertex_s *v);

    }
}

static vchain_s *vchain_s_ (void)
{
    vchain_s *c;
    
    c = calloc(1,sizeof(*c));
    if (!c)
        return NULL;
    c->mem = ~c->mem;
    return c;
}

int chain_insert (vchain_s *chunk, vertex_s *v)
{
    uint8_t i;
    uint8_t it;
    vchain_s *iter;
    
    if (!chunk) {
        chunk = vchain_s_();
        if (!chunk)
            return 0;
    }
    for (iter = chunk; iter->next; iter = iter->next);
    i = ffs(iter->mem);
    if (i) {
        iter->mem &= ~(1 << i);
        i--;
    }
    else {
        iter->next = vchain_s_();
        iter = iter->next;
        if (!iter)
            return 0;
        iter->mem = 0xef;
    }
    iter->chunk[i].v = v;
    return 1;
}