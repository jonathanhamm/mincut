#include "gparse.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>  

vhash_s vhash_;
static gtoken_s *stream_;

/*reads file into a buffer*/
static unsigned char *read_gfile(const char *fname);
static void tablegen (void);

/*tokenizing routines for graph data*/
static gtoken_s *gtoken_s_ (gtoken_s *node, unsigned char *lexeme, unsigned short type);

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
static void printbyte (uint8_t b);

int checkduplicates (wgraph_s *g, vertex_s *v1, vertex_s *v2)
{
  int i, j;
  
  for (i = 0; i < g->nvert; i++) {
    if (g->vtable[i] == v1) {
      for (j = 0; j < v1->nedges; j++) {
        if (v1->edges[j]->v1 == v1) {
          if (v1->edges[j]->v2 == v2)
            return 1;
        }
        else if (v1->edges[j]->v1 == v2){
          if (v1->edges[j]->v1 == v1)
            return 1;
        }
      }
    }
  }
  return 0;
}

wgraph_s *gparse (const unsigned char *file)
{
  wgraph_s *g;
  unsigned char *buf;
  
  buf = read_gfile (file);
  if (!buf)
    return NULL;
  if (!lex_ (buf))
    return NULL;
  g = parse_();
  if (!g)
    return NULL;
  return g;
}

unsigned char *read_gfile (const char *fname)
{
  FILE *f;
  size_t offset,
  bsize;
  unsigned char *buf;
  
  f = fopen(fname,"r");
  if(!f) {
    perror("Error Opening File");
    exit(EXIT_FAILURE);
  }
  buf = malloc(INITBUFSIZE);
  if (!buf)
    goto err_;
  for (bsize = INITBUFSIZE, offset = 0; (buf[offset] = (unsigned char)fgetc(f)) != UEOF; offset++) {
    if (offset == bsize-1) {
      bsize += INITBUFSIZE;
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

int vhashinsert (vertex_s *v, uint16_t index)
{
  uint16_t i;
  vrec_s *ptr;

  i = (unsigned long)v % VHTABLESIZE;
  if (vhash_.table[i].isoccupied) {
    ptr = malloc (sizeof(*ptr));
    if (!ptr)
      return -1;
    ptr->v = v;
    ptr->index = index;
    if (vhash_.table[i].isoccupied == 1) {
      ptr->next = NULL;
      vhash_.table[i].next = ptr;
    }
    else {
      ptr->next = vhash_.table[i].next;
      vhash_.table[i].next = ptr;
    }
  }
  else {
    vhash_.table[i].v = v;
    vhash_.table[i].index = index;
    vhash_.table[i].isoccupied = 1;
  }
  return 1;
}

uint16_t vgetindex (vertex_s *v)
{
  uint16_t i;
  vrec_s *ptr;
  
  i = (unsigned long)v % VHTABLESIZE;
  if (vhash_.table[i].v == v)
    return vhash_.table[i].index;
  else {
    for (ptr = vhash_.table[i].next; ptr; ptr = ptr->next) {
      if (ptr->v == v)
        return ptr->index;
    }
  }
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
  
  for (curr = NULL, bckptr = buf; *buf != UEOF;) {
    switch (*buf) {
      case ',':
        curr = gtoken_s_ (curr, ",", T_COMMA);
        buf++;
        break;
      case '=':
        curr = gtoken_s_ (curr, "=", T_EQU);
        buf++;
        break;
      case '{':
        curr = gtoken_s_ (curr, "{", T_OPENBRACE);
        buf++;
        break;
      case '}':
        buf++;
        curr = gtoken_s_ (curr, "}", T_CLOSEBRACE);
        break;
      default:
        if (*buf <= ' ')
          while(*++buf <= ' ');
        else if ((*buf >= 'A' && *buf <= 'Z') || (*buf >= 'a' && *buf <= 'z')) {
          for (bckptr = buf, ++buf; (*buf >= 'A' && *buf <= 'Z') || (*buf >= 'a' && *buf <= 'z')
               || (*buf >= '0' && *buf <= '9'); buf++) {
            if (buf - bckptr == MAXLEXLEN) {
              printf("Too Long ID: %15s", bckptr);
              goto err_;
            }
          }
          backup = *buf;
          *buf = '\0';
          curr = gtoken_s_ (curr, bckptr, T_ID);
          *buf = backup;
        } else if (*buf >= '0' && *buf <= '9') {
          for (bckptr = buf, buf++; (*buf >= '0' && *buf <= '9'); buf++) {
            if (buf - bckptr == MAXLEXLEN) {
              printf("Too Long ID: %15s", bckptr);
              goto err_;
            }
          }
          if (*buf == '.') {
            for (buf++; (*buf >='0' && *buf <= '9'); buf++) {
              if (buf - bckptr == MAXLEXLEN) {
                printf("Too Long ID: %15s", bckptr);
                goto err_;
              }
            }
          }
          backup = *buf;
          *buf = '\0';
          curr = gtoken_s_ (curr, bckptr, T_NUM);
          *buf = backup;
        } else {
          printf("Symbol Error %s %c\n", bckptr, *buf);
          goto err_;
        }
        break;
    }
  }
  curr = gtoken_s_ (curr, "$", T_EOF);
  return stream_;
err_:
  freetokens (stream_);
  stream_ = NULL;
  return NULL;
}

void freetokens (gtoken_s *list)
{
  gtoken_s *backup;
  
  while (list) {
    backup = list;
    list = list->next;
    free (backup);
  }
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


/*  Graph Grammar:
 *  <graph> => <nodelist> <edgelist> EOF
 *  <nodelist> => V={v <nodeparam>}
 *  <nodeparam> => ,v <nodeparam> | E
 *  <edgelist> => E={<e> <edgeparam> }
 *  <edgeparam> => ,<e> <edgeparam> | E
 *  <e> => {n,n,real}
 */
wgraph_s *parse_ (void)
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
  if (__GTNEXT()->type == T_EOF)
    printf("Parse Success!\n");
}

void pnodelist_ (wgraph_s *g)
{
  vertex_s *v;
  
  if (*(uint16_t *)stream_->lexeme == *(uint16_t *)"V") /*if(!strcmp(stream_->lexeme,"v"))*/
  if (__GTNEXT()->type == T_EQU)
  if (__GTNEXT()->type == T_OPENBRACE)
  if (__GTNEXT()->type == T_ID) {
    v = vertex_s_(stream_);
    vhashinsert (v, g->nvert);
    insert_vertex (g, v);
    pnodeparam_(g);
  }
  if (stream_->type == T_CLOSEBRACE)
    return;
}

void pnodeparam_ (wgraph_s *g)
{
  vertex_s *v;
  if (__GTNEXT()->type == T_COMMA)
  if (__GTNEXT()->type == T_ID) {
    v = vertex_s_(stream_);
    vhashinsert (v, g->nvert);
    insert_vertex (g, v);
    pnodeparam_(g);
  }
}

void pedgelist_ (wgraph_s *g)
{
  if (*(uint16_t *)__GTNEXT()->lexeme == *(uint16_t *)"E")   /*if (!strcmp(__GTNEXT()->lexeme, "E"))*/
  if (__GTNEXT()->type == T_EQU)
  if (__GTNEXT()->type == T_OPENBRACE)
    e_(g);
  pedgeparam_(g);
  if (stream_->type == T_CLOSEBRACE)
    return;
}

void pedgeparam_ (wgraph_s *g)
{
  while (__GTNEXT()->type == T_COMMA)
    e_(g);
}

void e_ (wgraph_s *g)
{
  float      weight;
  vertex_s    *v1,
              *v2;
  
  if (__GTNEXT()->type == T_OPENBRACE)
  if (__GTNEXT()->type == T_ID) {
    v1 = v_lookup (g, stream_->lexeme);
    if (!v1)
      return;
    if (__GTNEXT()->type == T_COMMA)
    if (__GTNEXT()->type == T_ID) {
      v2 = v_lookup (g, stream_->lexeme);
      if (!v2)
        return;
      if (__GTNEXT()->type == T_COMMA)
      if (__GTNEXT()->type == T_NUM) {
        weight = atof(stream_->lexeme);
        if (__GTNEXT()->type == T_CLOSEBRACE) {
          g->nedges++;
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
  
  //printf("Adding edge: %s %s %f\n",v1->name,v2->name, weight);
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
  static uint16_t cvtablesize;
  vertex_s **vtable;
  /*uint8_t index;
  vrecord_s *rec;
  
  rec = &graph->vtable[*(uint64_t *)v->name % _VTABLE_SIZE];
  if (!rec->isoccupied) {
    rec->v = v;
    rec->isoccupied = 1;
  } else if (rec->isoccupied >= 1)
    return chain_insert (&rec->chain, v);*/

  vtable = graph->vtable;
  if (!vtable) {
    vtable = calloc(INITTSIZE, sizeof(*vtable));
    if (!vtable)
      return -1;
    cvtablesize = INITTSIZE;
  } else {
    if (graph->nvert == cvtablesize) {
      cvtablesize += INITTSIZE;
      vtable = realloc (vtable, cvtablesize * sizeof(*vtable));
      if (!vtable)
        return -1;
    }
  }
  vtable[graph->nvert] = v;
  graph->vtable = vtable;
  graph->nvert++;
}

vertex_s *v_lookup (wgraph_s *graph, unsigned char *key)
{
  uint16_t i;
  uint8_t it;
  /*
  rec = &graph->vtable[*(uint64_t *)key % _VTABLE_SIZE];
  if (!__IDCMP(key,rec->v->name)) {
    if (rec->isoccupied && rec->isoccupied != 1) {
      for (iterator = rec->chain; iterator; iterator = iterator->next) {
        for (i = 0, it = ~iterator->mem; it; it &= ~(1 << i), i = ffs(it)-1) {
          if (__IDCMP(key,iterator->chunk[i].v->name))
            return iterator->chunk[i].v;
        }
      }
      if (!iterator)
        return NULL;
    }
    return NULL;
  }*/
  for (i = 0; i < graph->nvert; i++) {
    if (!strcmp(key,graph->vtable[i]->name))
      return graph->vtable[i];
  }
  return NULL;
}

void printbyte (uint8_t b)
{
  uint8_t i;
  
  for (i = 0; i < 8; i++)
    printf("%d",((b << i) & 0x80)>>7);
  printf("\n");
}