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

#define __IDCPY(DST,SRC)    *(uint64_t *)DST = *(uint64_t *)SRC; \
                            *(((uint64_t *)DST) + 1) = *(((uint64_t *)&SRC[0]) + 1)

#define __IDCMP(ID1,ID2)    (*(uint64_t *)ID1 == *(uint64_t *)ID2 && \
                            *(((uint64_t *)ID1) + 1) == *(((uint64_t *)ID2) + 1))

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
    edge_s **edges;
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
static wgraph_s *parse_ (void);
static void pgraph_ (wgraph_s *g);
static void pnodelist_ (wgraph_s *g);
static void pnodeparam_ (wgraph_s *g);
static void pedgelist_ (wgraph_s *g);
static void pedgeparam_ (wgraph_s *g);
static void e_ (wgraph_s *g);

/*graph data structure routines*/
static inline wgraph_s *wgraph_s_ (void);
static vertex_s *vertex_s_ (gtoken_s *tok);
static int addedge (vertex_s *v, edge_s *e);
static edge_s *edge_s_ (vertex_s *v1, vertex_s *v2, double weight);
static int insert_vertex (wgraph_s *graph, vertex_s *v);
static vertex_s *v_lookup (wgraph_s *graph, unsigned char *key);
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
static wgraph_s *parse_ (void)
{
    wgraph_s *g;
    
    g = wgraph_s_();
    if (!g)
        return NULL;
    pgraph_(g);
    return g;
}

void pgraph_ (wgraph_s *g)
{
    pnodelist_(g);
    pedgelist_(g);
    if (__GTNEXT()->type == _EOF)
        printf("Parse Success!\n");
}

void pnodelist_ (wgraph_s *g)
{
    if (!strcmp(stream_->lexeme, "V"))
    if (__GTNEXT()->type == _EQU)
    if (__GTNEXT()->type == _OPENBRACE)
    if (__GTNEXT()->type == _ID) {
        insert_vertex (g, vertex_s_(stream_));
        pnodeparam_(g);
    }
    if (stream_->type == _CLOSEBRACE)
        return;
}

void pnodeparam_ (wgraph_s *g)
{
    if (__GTNEXT()->type == _COMMA)
    if (__GTNEXT()->type == _ID) {
        insert_vertex (g, vertex_s_(stream_));
        pnodeparam_(g);
    }
}

void pedgelist_ (wgraph_s *g)
{
    if (!strcmp(__GTNEXT()->lexeme, "E"))
    if (__GTNEXT()->type == _EQU)
    if (__GTNEXT()->type == _OPENBRACE)
        e_(g);
    pedgeparam_(g);
    if (stream_->type == _CLOSEBRACE)
        return;
}

void pedgeparam_ (wgraph_s *g)
{
    if (__GTNEXT()->type == _COMMA) {
        e_(g);
        pedgeparam_(g);
    }
}

void e_ (wgraph_s *g)
{
    double      weight;
    vertex_s    *v1,
                *v2;
    
    if (__GTNEXT()->type == _OPENBRACE)
    if (__GTNEXT()->type == _ID) {
        v1 = v_lookup (g, stream_->lexeme);
        if (!v1)
            return;
        if (__GTNEXT()->type == _COMMA)
        if (__GTNEXT()->type == _ID) {
            v2 = v_lookup (g, stream_->lexeme);
            if (!v2)
                return;
            if (__GTNEXT()->type == _COMMA)
            if (__GTNEXT()->type == _NUM) {
                weight = atof(stream_->lexeme);
                if (__GTNEXT()->type == _CLOSEBRACE) {
                    edge_s_ (v1, v2, weight);
                    return;
                }
            }
        }
    }
}

/*graph data structure routines*/

inline wgraph_s *wgraph_s_ (void)
{
    return calloc(1,sizeof(wgraph_s));
}

vertex_s *vertex_s_ (gtoken_s *tok)
{
    vertex_s *v;
    
    v = calloc(1, sizeof(*v));
    if (!v)
        return NULL;
    __IDCPY (v->name, tok->lexeme);
    return v;
}

int addedge (vertex_s *v, edge_s *e)
{
    if (v->nedges)
        v->edges = realloc(v->edges, (v->nedges + 1) * sizeof(*v->edges));
    else
        v->edges = malloc(sizeof(*v->edges));
    if (!v->edges)
        return 0;
    v->edges[v->nedges] = e;
    v->nedges++;
    return 1;
}

edge_s *edge_s_ (vertex_s *v1, vertex_s *v2, double weight)
{
    edge_s *edge;
    
    edge = malloc(sizeof(*edge));
    if (!edge)
        return NULL;
    if (!(addedge(v1,edge) && addedge(v2,edge)))
        return NULL;
    edge->v1 = v1;
    edge->v2 = v2;
    edge->weight = weight;
    return edge;
}

int insert_vertex (wgraph_s *graph, vertex_s *v)
{
    uint8_t index;
    vrecord_s *rec;
    
    rec = &graph->vtable[*(uint64_t *)v->name % _VTABLE_SIZE];
    if (!rec->isoccupied) {
        rec->v = v;
        rec->isoccupied = 1;
    } else if (rec->isoccupied >= 1)
        return chain_insert (rec->chain, v);
}

vertex_s *v_lookup (wgraph_s *graph, unsigned char *key)
{
    uint8_t i;
    uint8_t it;
    vrecord_s *rec;
    vchain_s *iterator;
    
    rec = &graph->vtable[*(uint64_t *)key % _VTABLE_SIZE];
    if (!__IDCMP(key,rec->v->name)) {
        if (rec->isoccupied && rec->isoccupied != 1) {
            for (iterator = rec->chain; iterator; iterator = iterator->next) {
                for (i = 0, it = iterator->mem; it; it &= ~(1 << i), i = ffs(it)-1) {
                    if (__IDCMP(key,iterator->chunk[i].v->name))
                        return iterator->chunk[i].v;
                }
            }
            if (!iterator)
                return NULL;
        }
        return NULL;
    }
    return rec->v;
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
    } else {
        iter->next = vchain_s_();
        iter = iter->next;
        if (!iter)
            return 0;
        iter->mem = 0x7f; /* 0b01111111 */
    }
    iter->chunk[i].v = v;
    return 1;
}