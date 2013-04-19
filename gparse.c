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
  if(!f)
    THROW_EXCEPTION();
  buf = malloc(INITBUFSIZE);
  if (!buf)
    THROW_EXCEPTION();
  for (bsize = INITBUFSIZE, offset = 0; (buf[offset] = (unsigned char)fgetc(f)) != UEOF; offset++) {
    if (offset == bsize-1) {
      bsize += INITBUFSIZE;
      buf = realloc (buf, bsize);
      if (!buf)
        THROW_EXCEPTION();
    }
  }
  /*truncate buffer to EOF*/
  if (offset < bsize)
    buf = realloc (buf, offset+1);
  return buf;
  
exception_:
  perror("File Read Error");
  fclose(f);
  exit(EXIT_FAILURE);
}

int vhashinsert (vertex_s *v, uint16_t index)
{
  uint16_t i;
  vrec_s *ptr;

  i = (unsigned long)v % VHTABLESIZE;
  if (vhash_.table[i].isoccupied) {
    ptr = malloc (sizeof(*ptr));
    if (!ptr) {
      perror("Malloc Error");
      exit(EXIT_FAILURE);
    }
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
  
  stream_ = NULL;
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
              THROW_EXCEPTION();
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
              THROW_EXCEPTION();
            }
          }
          if (*buf == '.') {
            for (buf++; (*buf >='0' && *buf <= '9'); buf++) {
              if (buf - bckptr == MAXLEXLEN) {
                printf("Too Long ID: %15s", bckptr);
                THROW_EXCEPTION();
              }
            }
          }
          backup = *buf;
          *buf = '\0';
          curr = gtoken_s_ (curr, bckptr, T_NUM);
          *buf = backup;
        } else {
          printf("Symbol Error %c\n", *bckptr);
          THROW_EXCEPTION();
        }
        break;
    }
  }
  curr = gtoken_s_ (curr, "$", T_EOF);
  return stream_;
  
exception_:
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
  if (!tok) {
    perror ("Malloc Error");
    exit(EXIT_FAILURE);
  }
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
  if (!g) {
    perror("MALLOC ERROR");
    exit(EXIT_FAILURE);
  }
  pgraph_(g);
  return g;
}

void pgraph_ (wgraph_s *g)
{
  pnodelist_(g);
  pedgelist_(g);
  if (GTNEXT()->type != T_EOF) {
    printf ("Syntax Error: %s", stream_->lexeme);
    exit(EXIT_FAILURE);
  }
  printf ("Parse Success.\n");
}

void pnodelist_ (wgraph_s *g)
{
  vertex_s *v;
  
  if (*(uint16_t *)stream_->lexeme == *(uint16_t *)"V") /*if(!strcmp(stream_->lexeme,"v"))*/
  if (GTNEXT()->type == T_EQU)
  if (GTNEXT()->type == T_OPENBRACE)
  if (GTNEXT()->type == T_ID) {
    v = vertex_s_(stream_);
    vhashinsert (v, g->nvert);
    insert_vertex (g, v);
    pnodeparam_(g);
  }
  if (stream_->type == T_CLOSEBRACE)
    return;
  printf ("Syntax Error: %s", stream_->lexeme);
  exit(EXIT_FAILURE);
}

void pnodeparam_ (wgraph_s *g)
{
  vertex_s *v;
  
  if (GTNEXT()->type == T_COMMA)
  if (GTNEXT()->type == T_ID) {
    v = vertex_s_(stream_);
    vhashinsert (v, g->nvert);
    insert_vertex (g, v);
    pnodeparam_(g);
    return;
  }
}

void pedgelist_ (wgraph_s *g)
{
  if (*(uint16_t *)GTNEXT()->lexeme == *(uint16_t *)"E")   /*if (!strcmp(__GTNEXT()->lexeme, "E"))*/
  if (GTNEXT()->type == T_EQU)
  if (GTNEXT()->type == T_OPENBRACE)
    e_(g);
  pedgeparam_(g);
  if (stream_->type == T_CLOSEBRACE)
    return;
  printf ("Syntax Error: %s\n", stream_->lexeme);
  exit(EXIT_FAILURE);
}

void pedgeparam_ (wgraph_s *g)
{
  while (GTNEXT()->type == T_COMMA)
    e_(g);
}

void e_ (wgraph_s *g)
{
  float      weight;
  vertex_s    *v1,
              *v2;
  
  if (GTNEXT()->type == T_OPENBRACE)
  if (GTNEXT()->type == T_ID) {
    v1 = v_lookup (g, stream_->lexeme);
    if (!v1)
      THROW_EXCEPTION();
    if (GTNEXT()->type == T_COMMA)
    if (GTNEXT()->type == T_ID) {
      v2 = v_lookup (g, stream_->lexeme);
      if (!v2)
        THROW_EXCEPTION();
      if (GTNEXT()->type == T_COMMA)
      if (GTNEXT()->type == T_NUM) {
        weight = atof(stream_->lexeme);
        if (GTNEXT()->type == T_CLOSEBRACE) {
          g->nedges++;
          edge_s_ (v1, v2, weight);
          return;
        }
      }
    }
  }
  
exception_:
  printf ("Syntax Error: %s\n", stream_->lexeme);
  exit(EXIT_FAILURE);
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
  if (!v) {
    perror ("Malloc Error");
    exit(EXIT_FAILURE);
  }
  IDCPY (v->name, tok->lexeme);
  return v;
}

int addedge (vertex_s *v, edge_s *e)
{
  if (v->nedges)
    v->edges = realloc(v->edges, (v->nedges + 1) * sizeof(*v->edges));
  else
    v->edges = malloc(sizeof(*v->edges));
  if (!v->edges) {
    perror("Malloc Error");
    exit(EXIT_FAILURE);
  }
  v->edges[v->nedges] = e;
  v->nedges++;
  return 1;
}

edge_s *edge_s_ (vertex_s *v1, vertex_s *v2, double weight)
{
  edge_s *edge;
  
  edge = malloc(sizeof(*edge));
  if (!edge)
    THROW_EXCEPTION();
  if (!(addedge(v1,edge) && addedge(v2,edge)))
    THROW_EXCEPTION();
  edge->v1 = v1;
  edge->v2 = v2;
  edge->weight = weight;
  return edge;

exception_:
  perror("Malloc Error");
  exit(EXIT_FAILURE);
}

int insert_vertex (wgraph_s *graph, vertex_s *v)
{
  static uint16_t cvtablesize;
  vertex_s **vtable;

  vtable = graph->vtable;
  if (!vtable) {
    vtable = calloc(INITTSIZE, sizeof(*vtable));
    if (!vtable)
      THROW_EXCEPTION();
    cvtablesize = INITTSIZE;
  } else {
    if (graph->nvert == cvtablesize) {
      cvtablesize += INITTSIZE;
      vtable = realloc (vtable, cvtablesize * sizeof(*vtable));
      if (!vtable)
        THROW_EXCEPTION();
    }
  }
  vtable[graph->nvert] = v;
  graph->vtable = vtable;
  graph->nvert++;
  return 1;

exception_:
  perror("Malloc Error");
  exit(EXIT_FAILURE);
}

vertex_s *v_lookup (wgraph_s *graph, unsigned char *key)
{
  uint16_t i;
  uint8_t it;

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