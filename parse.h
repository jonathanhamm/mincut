/*
 gparse.h
 Author: Jonathan Hamm
 
 Description:
 
 This library is exclusively for routines that read and parse data so the genetic algorithm
 can use it. This contains functions for file io, parsing graph data, and parsing command line
 data. Command line data are commands sent to the GA during runtime from the terminal. These
 are mostly for printing results and setting parameters while the algorithm is running. This
 also contains functions for creating the graph data structure used by the GA.
 */

#ifndef _GPARSE_H_
#define _GPARSE_H_

#include <stdint.h>

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

#define throw_exception() goto exception_

#define IDCPY(DST,SRC)    *(uint64_t *)DST = *(uint64_t *)SRC; \
*(((uint64_t *)DST) + 1) = *(((uint64_t *)SRC) + 1)

#define GTNEXT() (stream_ = stream_->next)
#define GTPREV() (stream_ = stream_->prev)


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
extern gtoken_s *stream_;

/* Graph Parsing Routines */
extern wgraph_s *gparse (const unsigned char *file);
extern gtoken_s *lex (unsigned char *buf);
extern void freetokens (gtoken_s *list);

/*graph data structure routines*/
extern inline wgraph_s *wgraph_s_ (void);
extern vertex_s *vertex_s_ (gtoken_s *tok);
extern int addedge (vertex_s *v, edge_s *e);
extern edge_s *edge_s_ (vertex_s *v1, vertex_s *v2, double weight);
extern int insert_vertex (wgraph_s *graph, vertex_s *v);
extern int vhashinsert (vertex_s *v, uint16_t index);
extern uint16_t vgetindex (vertex_s *v);
extern void printgraph (wgraph_s *g);
extern void printsastatus (void);

/* Command line parsing routines */
extern void cgeparse (void);
extern void csaparse (void);

#endif
