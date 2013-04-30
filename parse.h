/*
 parse.h
 Author: Jonathan Hamm
 
 Description:
 
 Library implemented by parse.c
 
 This is a library for all the parsing operations needed 
 for the genetic algorithm, simulated annealing, and foolish
 hill climbing. 
 */

#ifndef _GPARSE_H_
#define _GPARSE_H_

#include <stdint.h>

#define CHUNK_SIZE 8
#define UEOF (unsigned char)EOF
#define INITBUFSIZE 256
#define MAXLEXLEN 15
#define INITTSIZE 32
#define VHTABLESIZE 233

/*Token Types*/
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

typedef struct vhash_s vhash_s;
typedef struct vrec_s vrec_s;

/* Graph Data Structure */
struct wgraph_s
{
    uint32_t  nedges;
    uint16_t  nvert;
    vertex_s  **vtable;
};

/* Token structure */
struct gtoken_s
{
    unsigned short type;
    unsigned char lexeme[MAXLEXLEN + 1];
    gtoken_s *prev;
    gtoken_s *next;
};

/* Graph vertex Structure */
struct vertex_s
{
    unsigned char name[MAXLEXLEN + 1];
    uint32_t nedges;
    edge_s **edges;
};

/* Graph Edge Structure */
struct edge_s
{
    double weight;
    vertex_s *v1;
    vertex_s *v2;
};

/* Record in hash table for vertex indices */
struct vrec_s
{
    vertex_s *v;
    uint16_t index;
    union {
        unsigned long isoccupied;
        vrec_s *next;
    };
};

/* Hash Table for vertex indices. */
struct vhash_s
{
    vrec_s table[VHTABLESIZE];
};

/* Global Variables */
extern vhash_s  vhash_;     //vertex hash
extern gtoken_s *stream_;   //linked list of tokens (token stream)

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