#ifndef _GPARSE_H_
#define _GPARSE_H_

#include <stdint.h>

#define MAPNODE_SIZE 4
#define CHUNK_SIZE 8
#define UEOF (unsigned char)EOF
#define INITBUFSIZE 256
#define MAXLEXLEN 15
#define INITTSIZE 32

/*TTYPES*/
#define T_ID         1
#define T_NUM        2
#define T_OPENBRACE  3
#define T_CLOSEBRACE 4
#define T_COMMA      5
#define T_EQU        6
#define T_EOF        7

#define __IDCPY(DST,SRC)    *(uint64_t *)DST = *(uint64_t *)SRC; \
*(((uint64_t *)DST) + 1) = *(((uint64_t *)SRC) + 1)

#define __IDCMP(ID1,ID2)    (*(uint64_t *)ID1 == *(uint64_t *)ID2 && \
*(((uint64_t *)ID1) + 1) == *(((uint64_t *)ID2) + 1))

#define __GTNEXT() (stream_ = stream_->next)
#define __GTPREV() (stream_ = stream_->prev)


typedef struct gnode_s gnode_s;
typedef struct wgraph_s wgraph_s;
typedef struct gtoken_s gtoken_s;
typedef struct vertex_s vertex_s;
typedef struct edge_s edge_s;

#define VHTABLESIZE 233
typedef struct vhash_s vhash_s;
typedef struct vrec_s vrec_s;

struct wgraph_s
{
  uint16_t  nedges;
  uint16_t  nvert;
  vertex_s  **vtable;
};


struct gtoken_s
{
  unsigned short type;
  unsigned char lexeme[MAXLEXLEN + 1];
  gtoken_s *prev;
  gtoken_s *next;
};

struct vertex_s
{
  unsigned char name[MAXLEXLEN + 1];
  uint16_t nedges;
  edge_s **edges;
};

struct edge_s
{
  float weight;
  vertex_s *v1;
  vertex_s *v2;
};

struct vrec_s
{
  vertex_s *v;
  uint16_t index;
  union {
    unsigned long isoccupied;
    vrec_s *next;
  };
};

struct vhash_s
{
  vrec_s table[VHTABLESIZE];
};

extern vhash_s vhash_;

extern gtoken_s *lex_ (unsigned char *buf);
extern void freetokens (gtoken_s *list);

/*graph data structure routines*/
extern wgraph_s *gparse (const unsigned char *file);

extern int vhashinsert (vertex_s *v, uint16_t index);
extern uint16_t vgetindex (vertex_s *v);

extern void printgraph (wgraph_s *g);

#endif
