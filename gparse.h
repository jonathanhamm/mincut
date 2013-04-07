#ifndef _GPARSE_H_
#define _GPARSE_H_

#include <stdint.h>

#define _MAPNODE_SIZE 4
#define _CHUNK_SIZE 8
#define _UEOF (unsigned char)EOF
#define _INITBUFSIZE 256
#define _MAXLEXLEN 15
#define _INITTSIZE 32

/*TTYPES*/
#define _ID         1
#define _NUM        2
#define _OPENBRACE  3
#define _CLOSEBRACE 4
#define _COMMA      5
#define _EQU        6
#define _EOF        7

#define __IDCPY(DST,SRC)    *(uint64_t *)DST = *(uint64_t *)SRC; \
*(((uint64_t *)DST) + 1) = *(((uint64_t *)SRC) + 1)

#define __IDCMP(ID1,ID2)    (*(uint64_t *)ID1 == *(uint64_t *)ID2 && \
*(((uint64_t *)ID1) + 1) == *(((uint64_t *)ID2) + 1))

#define __GTNEXT() (stream_ = stream_->next)
#define __GTPREV() (stream_ = stream_->prev)


typedef struct gnode_s gnode_s;
typedef struct wgraph_s wgraph_s;
typedef struct gtoken_s gtoken_s;
typedef struct vrecord_s vrecord_s;
typedef struct vertex_s vertex_s;
typedef struct vchain_s vchain_s;
typedef struct edge_s edge_s;
typedef struct wmap_s wmap_s;

struct vrecord_s
{
  vertex_s *v;
  union {
    unsigned long isoccupied;
    vchain_s *chain;
  };
};

struct wgraph_s
{
  uint16_t  nedges;
  uint16_t  nvert;
  vertex_s  **vtable;
};


struct gtoken_s
{
  unsigned short type;
  unsigned char lexeme[_MAXLEXLEN + 1];
  gtoken_s *prev;
  gtoken_s *next;
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
  uint16_t nedges;
  edge_s **edges;
};

struct edge_s
{
  float weight;
  vertex_s *v1;
  vertex_s *v2;
};

struct wmap_s
{
  float map[_MAPNODE_SIZE][256];
  wmap_s *next;
};

/*graph data structure routines*/
extern wgraph_s *gparse (const unsigned char *file);
extern void printgraph (wgraph_s *g);

#endif
