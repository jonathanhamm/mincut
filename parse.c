#include "parse.h"
#include "bis.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>  

#define COM_ERR   0
#define COM_PROB  1
#define COM_OP    2
#define COM_X     3
#define COM_SEL   4
#define COM_K     5

vhash_s vhash_;
gtoken_s *stream_;
pool_s *pool_;


/* reads file into a buffer */
static unsigned char *read_gfile(const char *fname);
static void tablegen (void);

/* tokenizing routines for graph data */
static gtoken_s *gtoken_s_ (gtoken_s *node, unsigned char *lexeme, unsigned short type);

/* graph parsing routines */
static wgraph_s *parse_ (void);
static void pgraph_ (wgraph_s *g);
static void pnodelist_ (wgraph_s *g);
static void pnodeparam_ (wgraph_s *g);
static void pedgelist_ (wgraph_s *g);
static void pedgeparam_ (wgraph_s *g);
static void e_ (wgraph_s *g);

/* graph data structure routines */
static vertex_s *v_lookup (wgraph_s *graph, unsigned char *key);
static void printbyte (uint8_t b);

/* Commandline Parsing Routines */
static int  p_op (void);
static int  p_mutate (void);
static void p_show (void);
static int p_feasible (void);


wgraph_s *gparse (const unsigned char *file)
{
  wgraph_s *g;
  unsigned char *buf;
  
  buf = read_gfile (file);
  if (!buf)
    return NULL;
  if (!lex (buf))
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
      perror("Heap Allocation Error");
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
gtoken_s *lex (unsigned char *buf)
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
    perror ("Heap Allocation Error");
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
 *  <nodeparam> => ,v <nodeparam> | epsilon
 *  <edgelist> => E={<e> <edgeparam> }
 *  <edgeparam> => ,<e> <edgeparam> | epsilon
 *  <e> => {n,n,real}
 */
wgraph_s *parse_ (void)
{
  wgraph_s *g;
  
  g = wgraph_s_();
  if (!g) {
    perror("Heap Allocation Error");
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
    perror ("Heap Allocation Error");
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
    perror("Heap Allocation Error");
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
  perror("Heap Allocation Error");
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
  perror("Heap Allocation Error");
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

/*
 Commandline Grammar
 
 <cparse>=>
 status | exit | quit | q | Q
 |
 set <op> numreal | get <op> | show <show>
 
 <op>=>
 mutate <mutate> | cross | select | k
 
 <mutate>=>
 op | prob
 
 <selection>=>
 selection id
 
 <show>=>
 numint | best <feasible>
 
 <feasible>=> feasible | epsilon
 */
void cparse (void)
{
  int result, val;
  
  if (
      !strcmp (stream_->lexeme, "exit")  ||
      !strcmp (stream_->lexeme, "quit")  ||
      !strcmp (stream_->lexeme, "Quit")  ||
      !strcmp (stream_->lexeme, "q")     ||
      !strcmp (stream_->lexeme, "Q")
      )
  {
    printf("Final:\n");
    printstatus ();
    kill(getppid(), SIGQUIT);
    exit(EXIT_SUCCESS);
  }
  else if (!strcmp(stream_->lexeme, "set")) {
    GTNEXT();
    result = p_op ();
    if (stream_->type == T_NUM) {
      val = atoi (stream_->lexeme);
      GTNEXT();
      switch (result) {
        case COM_OP:
          if (val <= 1) {
            pool_->mutate = mutate1;
            printf("Mutate operator now set to: 1\n");
          }
          else {
            pool_->mutate = mutate2;
            printf("Mutate operator now set to: 2\n");
          }
          break;
        case COM_PROB:
          if (val < 0) val = 0;
          else if (val > 100) val = 100;
          pool_->mutateprob = (uint8_t)val;
          printf("Mutate probability now set to %d\n", val);
          break;
        case COM_X:
          if (val <= 1) {
            pool_->cross = uniform_cr;
            printf("Now using uniform crossover\n");
          }
          else {
            pool_->cross = npoint_cr;
            printf("Now using n-point crossover, n = %d\n", CR_N);
          }
          break;
        case COM_SEL:
          if (val <= 1) {
            pool_->select = roulette_sf;
            printf("Now using roulette selection.");
          }
          else if (val == 2) {
            pool_->select = rank_sf;
            printf("Now using rank selection\n");
          }
          else {
            pool_->select = tournament_sf;
            printf("Now using tournament selection\n");
          }
          break;
        case COM_K:
          pool_->k = (uint8_t)val;
          printf ("Now using k value (for tournament selction): %d.\n", val);
          break;
        default:
          break;
      }
    }
    else
      printf ("Command Line Error: Expected number, but got '%s'\n", stream_->lexeme);
  }
  else if (!strcmp(stream_->lexeme, "get")) {
    GTNEXT();
    result = p_op();
    switch (result) {
      case COM_OP:
        if (pool_->mutate == mutate1)
          printf("Currently using mutation operator 1.\n");
        else
          printf("Currently using mutation operator 2.\n");
        break;
      case COM_PROB:
        printf("Current mutation probability is %u%%\n", pool_->mutateprob);
        break;
      case COM_X:
        if (pool_->cross == uniform_cr)
          printf("Currently using uniform crossover.\n");
        else
          printf ("Currently using n-point crossover, n = %d\n", CR_N);
        break;
      case COM_SEL:
        if (pool_->select == roulette_sf)
          printf("Currently using roulette selection.\n");
        else if (pool_->select == rank_sf)
          printf("Currently using rank selection.\n");
        else
          printf("Currently using tournament selection.\n");
        break;
      case COM_K:
        printf("Current k value (for tournament selection): %d.\n", pool_->k);
        break;
      default:
        break;
    }
    
  }
  else if (!strcmp(stream_->lexeme, "show")) {
    GTNEXT();
    p_show();
  }
  else if (!strcmp (stream_->lexeme, "status")) {
    GTNEXT();
    printstatus();
  }
  else
    printf ("Command Line Error: Unrecognized: '%s'\n", stream_->lexeme);
}

int p_op (void)
{
  if (!strcmp(stream_->lexeme, "mutate")) {
    GTNEXT();
    return p_mutate();
  }
  else if (!strcmp(stream_->lexeme, "cross")) {
    GTNEXT();
    return COM_X;
  }
  else if (!strcmp(stream_->lexeme, "select")) {
    GTNEXT();
    return COM_SEL;
  }
  else if (!strcmp(stream_->lexeme, "k")) {
    GTNEXT();
    return COM_K;
  }
  printf ("Command Line Error: Expected 'mutate', 'cross', 'select', or 'k' but got '%s'.\n", stream_->lexeme);
  return COM_ERR;
}

int p_mutate (void)
{
  if (!strcmp (stream_->lexeme, "op")) {
    GTNEXT();
    return COM_OP;
  }
  else if (!strcmp(stream_->lexeme, "prob")) {
    GTNEXT();
    return COM_PROB;
  }
  printf ("Command Line Error: Expected 'op' or 'prob' but got: '%s'.\n", stream_->lexeme);
  return COM_ERR;
}

void p_show (void)
{
  int result, index;
  
  if (stream_->type == T_NUM) {
    index = atoi(stream_->lexeme);
    GTNEXT();
    if (index <= 0 || index > POOLSIZE)
      printf ("Value %d out of range. Range is 1 to %d.\n", index, POOLSIZE);
    else
      printsolution(--index);
  }
  else if (!strcmp (stream_->lexeme, "best")) {
    GTNEXT();
    result = p_feasible ();
    if (!result)
      printsolution (POOLSIZE-1);
    else if (result == 1) {
      for (index = POOLSIZE-1; index >= 0
           && pool_->rbuf[index].ptr != pool_->bestfeasible; index--);
      printsolution (index);
    }
    
  }
  else
    printf ("Expected number or 'best', but got '%s'.\n", stream_->lexeme);
}

int p_feasible (void)
{
  int index;
  
  if (!strcmp (stream_->lexeme, "feasible")) {
    GTNEXT();
    return 1;
  }
  else if (stream_->type != T_EOF) {
    printf ("Expected 'feasible' or nothing, but got '%s'.\n", stream_->lexeme);
    return -1;
  }
  return 0;
}