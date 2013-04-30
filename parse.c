/*
 parse.c
 Author: Jonathan Hamm
 
 Description: 
 
 This file only contains code for parsing. None of the actual 
 genetic algorithm, simulated annealing, or foolish hill climbing
 code is here. However, this file does have code for building and 
 defining the data structures used by those algorithms. 
 */

#include "parse.h"
#include "bis.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Genetic Algorithm Commandline Constants */
#define COM_ERR     -1
#define COM_PROB    1
#define COM_OP      2
#define COM_X       3
#define COM_SEL     4
#define COM_K       5

/* Simulated Annealing Commandline Constants */
#define COM_T       1
#define COM_ITER    2
#define COM_ALPHA   3
#define COM_BETA    4
#define COM_PERTURB 5

/* Global Variables */
vhash_s     vhash_;     //hash (declared extern in parse.h)
gtoken_s    *stream_;   //token stream
pool_s      *pool_;     //pool of chromosomes (declared extern in bis.h, and linked in bis.c)

/* reads file into a buffer */
static unsigned char *read_gfile (const char *fname);

/* Token "constructor" */
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

/* Genetic Algorithm Commandline Parsing Routines */
static int  p_op (void);
static int  p_mutate (void);
static void p_show (void);
static int p_feasible (void);

/* Simulated Annealing Comandline Parsing Routines */
static int saparam (void);
static int sashow (void);

/*
 Function invoked for parsing a graph file.
 
 @param file    Name of file containing graph data. 
 @return        Returns a graph data structure.
 */ 
wgraph_s *gparse (const unsigned char *file)
{
    wgraph_s        *g;
    unsigned char   *buf;
    
    buf = read_gfile (file);
    if (!buf)
        return NULL;
    if (!lex (buf)) {
        perror("Parsing Graph Failed");
        exit(EXIT_FAILURE);
    }
    free(buf);
    g = parse_();
    if (!g)
        return NULL;
    return g;
}

/* 
 Reads a file into a char buffer. 
 
 @param fname   Name of file to read. 
 @return        Returns a pointer to the buffer holding the file data.
 */
unsigned char *read_gfile (const char *fname)
{
    FILE            *f;
    size_t          offset, bsize;
    unsigned char   *buf;
    
    f = fopen(fname,"r");
    if(!f)
        throw_exception();
    buf = malloc(INITBUFSIZE);
    if (!buf)
        throw_exception();
    for (bsize = INITBUFSIZE, offset = 0; (buf[offset] = (unsigned char)fgetc(f)) != UEOF; offset++) {
        if (offset == bsize-1) {
            bsize *= 2;
            buf = realloc (buf, bsize);
            if (!buf)
                throw_exception();
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

/*
 Inserts an index into a hash table. The index 
 corresponds to the index of a vertex in the 
 graph structure. The index is hashed on its 
 vertex's structure pointer for lookup. The 
 pointer/hash key is located when looking
 through the first set's edges. 
 
 @param v       The vertex structur pointer which is the
                hash key. 
 @param index   The index of the vertex. 
 @return        Returns 1 on success. 
 */
int vhashinsert (vertex_s *v, uint16_t index)
{
    uint16_t    i;
    vrec_s      *ptr;
    
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

/* 
 Looks at hash table to find the index of a 
 vertex in the graph structure. 
 
 @param v   Pointer to vertex structure that 
            is the lookup key.
 @return    Returns the index of the vertex in
            the graph structure. 
 */
uint16_t vgetindex (vertex_s *v)
{
    uint16_t    i;
    vrec_s      *ptr;
    
    i = (unsigned long)v % VHTABLESIZE;
    if (vhash_.table[i].v == v)
        return vhash_.table[i].index;
    else {
        for (ptr = vhash_.table[i].next; ptr; ptr = ptr->next) {
            if (ptr->v == v)
                return ptr->index;
        }
    }
    return i;
}

/*
 Small lexical analyzer used for tokenizing 
 a buffer of a data. This is used for tokenizing
 graph data from the input file and tokenizing
 commands entered by the user at runtime. 
 
 Lexer tokenizes based on the following regex: 
    token:  id | num
    id:     (a...Z)+ (a...Z | 0...9)*
    num:    magnit | - magnit
    magnit: (0...9)+ (dot (0...9)*)? | (dot (0...9)+)
 
 @param     buf     Pointer to the buffer that is tokenized.
 @return            Returns a pointer to a linked list of 
                    tokens on success, and NULL if there 
                    is a lexical error.
 */
gtoken_s *lex (unsigned char *buf)
{
    unsigned char   backup, gotnum,
                    *bckptr;
    gtoken_s        *curr;
    
    stream_ = NULL;
    gotnum = 0;
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
                            throw_exception();
                        }
                    }
                    backup = *buf;
                    *buf = '\0';
                    curr = gtoken_s_ (curr, bckptr, T_ID);
                    *buf = backup;
                } else if ((*buf >= '0' && *buf <= '9') || *buf == '.' || *buf == '-') {
                    bckptr = buf;
                    if (*buf >= '0' && *buf <= '9')
                        gotnum = 1;
                    else if (*buf == '-') {
                        buf++;
                        if ((*buf < '0' || *buf > '9') && *buf != '.') {
                            printf ("Symbol Error %.15s\n", bckptr);
                            throw_exception();
                        }
                        else if ((*buf >= '0' && *buf <= '9'))
                            gotnum = 1;
                    }
                    if (*buf != '.') {
                        for (buf++; (*buf >= '0' && *buf <= '9'); buf++) {
                            gotnum = 1;
                            if (buf - bckptr == MAXLEXLEN) {
                                printf("Too Long ID: %.15s", bckptr);
                                throw_exception();
                            }
                        }
                    }
                    if (*buf == '.') {
                        for (buf++; (*buf >='0' && *buf <= '9'); buf++) {
                            gotnum = 1;
                            if (buf - bckptr == MAXLEXLEN) {
                                printf("Too Long ID: %.15s", bckptr);
                                throw_exception();
                            }
                        }
                    }
                    if (gotnum) {
                        backup = *buf;
                        *buf = '\0';
                        gotnum = 0;
                        curr = gtoken_s_ (curr, bckptr, T_NUM);
                        *buf = backup;
                    }
                    else {
                        printf("Symbol Error %c\n", *bckptr);
                        throw_exception();
                    }
                } else {
                    printf("Symbol Error %c\n", *bckptr);
                    throw_exception();
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

/*
 Frees tokens (created by lex) from memory. 
 
 @param list    List of tokens to free. 
 */
void freetokens (gtoken_s *list)
{
    gtoken_s *backup;
    
    while (list) {
        backup = list;
        list = list->next;
        free (backup);
    }
}

/*
 "Constructs" a token and adds it to the list of nodes used
 by a parser.
 
 @param node    The tail of the list of tokens the new token
                will be appended to.
 @param lexeme  Pointer to the lexeme that the new token will 
                contain.
 @param type    The type of the token.
 @return        Returns a pointer to the new token "constructed".
 */
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

/*
 Parser for graph data. Parser adheres to
 the following grammar:
 
 Graph Parsing Grammar:
 <graph>=> 
    <nodelist> <edgelist> EOF
 
 <nodelist>=> 
    V = {id <nodeparam>}
 
 <nodeparam>=> 
    , id <nodeparam> | epsilon
 
 <edgelist>=> 
    E = {<e> <edgeparam> }
 
 <edgeparam>=> 
    , <e> <edgeparam> | epsilon
 
 <e>=> 
    {id, id, num}
 
 @return    Returns a pointer to the graph
            data structure built by the parser.
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
        printf ("Syntax Error: expected nothing, but got: %s", stream_->lexeme);
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
    printf ("Syntax Error while parsing nodes: %s", stream_->lexeme);
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
    printf ("Syntax Error while parsing in edge list (edge #%d): %s\n", g->nedges, stream_->lexeme);
    exit(EXIT_FAILURE);
}

void pedgeparam_ (wgraph_s *g)
{
    while (GTNEXT()->type == T_COMMA)
        e_(g);
}

void e_ (wgraph_s *g)
{
    double      weight;
    vertex_s    *v1,
                *v2;
    
    if (GTNEXT()->type == T_OPENBRACE)
    if (GTNEXT()->type == T_ID) {
        v1 = v_lookup (g, stream_->lexeme);
        if (!v1) {
            printf ("Edge: %s does not exist in node set.\n", stream_->lexeme);
            exit(EXIT_FAILURE);
        }
        if (GTNEXT()->type == T_COMMA)
        if (GTNEXT()->type == T_ID) {
            v2 = v_lookup (g, stream_->lexeme);
            if (!v2) {
                printf ("Edge: %s does not exist in node set.\n", stream_->lexeme);
                exit(EXIT_FAILURE);
            }
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
    printf ("Syntax Error while parsing an edge (edge #%d): %s\n", g->nedges, stream_->lexeme);
    exit(EXIT_FAILURE);
}

/*
 "Constructor" for graph structure. 
 
 @return    Returns a zeroed block of memory. 
 */
inline wgraph_s *wgraph_s_ (void)
{
    return calloc(1,sizeof(wgraph_s));
}

/*
 Vertex "constructor". Constructs a vertex for 
 the graph using a token. 
 
 @param tok Pointer to the token used for the vertex. 
 @return    Returns a pointer to the vertex structure.
 */
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

/*
 Adds an edge to a Vertex. 
 
 @param v   Vertex to attach edge to. 
 @param e   Edge to attach to vertex. 
 @return    Returns 1 on success
 */
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

/*
 Edge "constructor". Constructs an edge from 2 vertices with 
 a given weight. 
 
 @param v1      First vertex this edge connects. 
 @param v2      Second vertex this edge connects.
 @param weight  Weight of the edge.
 */
edge_s *edge_s_ (vertex_s *v1, vertex_s *v2, double weight)
{
    edge_s *edge;
    
    edge = malloc(sizeof(*edge));
    if (!edge)
        throw_exception();
    if (!(addedge(v1,edge) && addedge(v2,edge)))
        throw_exception();
    edge->v1 = v1;
    edge->v2 = v2;
    edge->weight = weight;
    return edge;
    
exception_:
    perror("Heap Allocation Error");
    exit(EXIT_FAILURE);
}

/*
 Inserts a vertex into the graph data strucure. 
 
 @param graph   The graph structure to insert the 
                vertex into.
 @param v       The vertex being inserted into the  
                graph.
 @return        Returns 1 on success.
 */
int insert_vertex (wgraph_s *graph, vertex_s *v)
{
    static uint16_t cvtablesize;
    vertex_s        **vtable;
    
    vtable = graph->vtable;
    if (!vtable) {
        vtable = calloc(INITTSIZE, sizeof(*vtable));
        if (!vtable)
            throw_exception();
        cvtablesize = INITTSIZE;
    } else {
        if (graph->nvert == cvtablesize) {
            cvtablesize += INITTSIZE;
            vtable = realloc (vtable, cvtablesize * sizeof(*vtable));
            if (!vtable)
                throw_exception();
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

/*
 Searches for a vertex in the vertex table of the graph. 
 The search key is the vertex's lexeme. 
 
 @param graph   Graph to search. 
 @param key     The search key, which is the vertex's lexeme. 
 @return        Returns a pointer to the vertex if found, 
                otherwise returns NULL if not found.
 */
vertex_s *v_lookup (wgraph_s *graph, unsigned char *key)
{
    uint16_t    i;
    
    for (i = 0; i < graph->nvert; i++) {
        if (!strcmp(key, graph->vtable[i]->name))
            return graph->vtable[i];
    }
    return NULL;
}

/*
 Parser for the genetic algorithm's commands 
 entered at runtime. Parser adheres to the 
 following grammar:

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
 
 <feasible>=> 
    feasible | epsilon
 */
void cgeparse (void)
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
        if (GTNEXT()->type != T_EOF) {
            printf ("Expected nothing, but got '%s'.\n", stream_->lexeme);
            return;
        }
        printf("Final:\n");
        printgestatus ();
        kill(getppid(), SIGQUIT);
        exit(EXIT_SUCCESS);
    }
    else if (!strcmp(stream_->lexeme, "set")) {
        GTNEXT();
        result = p_op();
        if (stream_->type == T_NUM) {
            val = atoi (stream_->lexeme);
            GTNEXT();
            if (stream_->type != T_EOF) {
                printf ("Expected nothing, but got '%s'.\n", stream_->lexeme);
                return;
            }
            switch (result) {
                case COM_OP:
                    if (val <= 1) {
                        pool_->mutate = mutate1;
                        printf("Mutate operator now set to: 1\n");
                    }
                    else if (val == 2){
                        pool_->mutate = mutate2;
                        printf("Mutate operator now set to: 2\n");
                    }
                    else {
                        pool_->mutate = pairwise_ex;
                        printf("Mutate operator now set to: 3 (pairwise exchange perturbation function)\n");
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
                        printf("Now using roulette selection.\n");
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
        if (stream_->type != T_EOF) {
            printf ("Expected nothing, but got '%s'.\n", stream_->lexeme);
            return;
        }
        switch (result) {
            case COM_OP:
                if (pool_->mutate == mutate1)
                    printf("Currently using mutation operator 1.\n");
                else if (pool_->mutate == mutate2)
                    printf("Currently using mutation operator 2.\n");
                else
                    printf("Currently using mutation operator 3 (pairwise exchange perturbation function).\n");
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
        printgestatus();
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
            printsolution(--index, NULL);
    }
    else if (!strcmp (stream_->lexeme, "best")) {
        GTNEXT();
        result = p_feasible();
        if (stream_->type != T_EOF) {
            printf ("Expected nothing, but got '%s'.\n", stream_->lexeme);
            return;
        }
        if (!result)
            printsolution (POOLSIZE-1, NULL);
        else if (result == 1) {
            for (index = POOLSIZE-1; index >= 0
                 && pool_->rbuf[index].ptr != pool_->bestfeasible; index--);
            printsolution (index, NULL);
        }
    }
    else
        printf ("Expected number or 'best', but got '%s'.\n", stream_->lexeme);
}

int p_feasible (void)
{
    if (!strcmp (stream_->lexeme, "feasible")) {
        GTNEXT();
        return 1;
    }
    else if (stream_->type != T_EOF) {
        printf ("Expected 'feasible' or nothing, but got '%s'.\n", stream_->lexeme);
        return COM_ERR;
    }
    return 0;
}

/*
 Parser for commands entered at runtime for simulated
 annealing and foolish hill climbing. The parser adheres
 to the following grammar:
 
 <cparse>=>
    status | exit | quit | q | Q
    |
    set <saparam> num
    |
    get <saparam>
    |
    show <show>
 
 <saparam>=>
    temp | alpha | iter | beta | perturb
 
 <show>=>
    best <feasible>
    |
    epsilon
 
 <feasible>=> 
    feasible | epsilon
 */
void csaparse (void)
{
    int     result;
    double  val;
    
    if (
        !strcmp (stream_->lexeme, "exit")  ||
        !strcmp (stream_->lexeme, "quit")  ||
        !strcmp (stream_->lexeme, "Quit")  ||
        !strcmp (stream_->lexeme, "q")     ||
        !strcmp (stream_->lexeme, "Q")
        )
    {
        if (GTNEXT()->type != T_EOF) {
            printf ("Expected nothing, but got '%s'.\n", stream_->lexeme);
            return;
        }
        printf("Final:\n");
        printsastatus ();
        kill(getppid(), SIGQUIT);
        exit(EXIT_SUCCESS);
    }
    else if (!strcmp(stream_->lexeme, "status")) {
        GTNEXT();
        printsastatus();
    }
    else if (!strcmp(stream_->lexeme, "set")) {
        result = saparam();
        GTNEXT();
        if (stream_->type == T_NUM && stream_->next->type)
            val = atof (stream_->lexeme);
        else {
            printf ("Expected number, but got '%s'.\n", stream_->lexeme);
            return;
        }
        if (GTNEXT()->type != T_EOF) {
            printf ("Expected nothing, but got '%s'.\n", stream_->lexeme);
            return;
        }
        switch (result) {
            case COM_T:
                pool_->T = val;
                printf("Temperature now set to: %f\n", val);
                break;
            case COM_ITER:
                pool_->iterations = val;
                printf("Number of iterationts now set to: %f\n", val);
                break;
            case COM_ALPHA:
                pool_->alpha = val;
                printf("Alpha now set to: %f\n", val);
                break;
            case COM_BETA:
                pool_->beta = val;
                printf("Beta now set to: %f\n", val);
                break;
            case COM_PERTURB:
                if (val <=1) {
                    pool_->perturb = mutate1;
                    printf("Set Perturbation Function to mutate function 1.\n");
                }
                else if (val == 2) {
                    pool_->perturb = mutate2;
                    printf("Set Perturbation Function to mutate function 2.\n");
                }
                else {
                    pool_->perturb = pairwise_ex;
                    printf("Set Perturbation Function to mutate function 3 (pairwise exchange).\n");
                }
                break;
            default:
                break;
        }
    }
    else if (!strcmp(stream_->lexeme, "get")) {
        result = saparam();
        if (GTNEXT()->type != T_EOF) {
            printf ("Expected nothing, but got '%s'.\n", stream_->lexeme);
            return;
        }
        switch (result) {
            case COM_T:
                printf ("Temperature currently set to: %f\n", pool_->T);
                break;
            case COM_ITER:
                printf ("Number of iterations currently set to: %f\n", pool_->iterations);
                break;
            case COM_ALPHA:
                printf ("Alpha currently set to: %f\n", pool_->alpha);
                break;
            case COM_BETA:
                printf ("Beta currently set to: %f\n", pool_->beta);
                break;
            case COM_PERTURB:
                if (pool_->perturb == mutate1)
                    printf("Current Perturbation Function is: mutate function 1.\n");
                else if (pool_->perturb == mutate1)
                    printf("Current Perturbation Function is: mutate function 2.\n");
                else 
                    printf("Current Perturbation Function is: mutate function 3 (pairwise exchange).\n");
                break;
            default:
                break;
        }
    }
    else if (!strcmp (stream_->lexeme, "show")) {
        GTNEXT();
        result = sashow();
        if (result == COM_ERR)
            return;
        if (stream_->type != T_EOF) {
            printf ("Expected nothing, but got '%s'.\n", stream_->lexeme);
            return;
        }
        switch (result) {
            case 0:
                printsolution (SIMA_curr, NULL);
                break;
            case 1:
                printsolution (SIMA_curr, NULL);
                break;
            case 2:
                printsolution (-1, pool_->bestfeasible);
            default:
                break;
        }
    }
    else
        printf ("Command Line Error: Unrecognized: '%s'\n", stream_->lexeme);
}

int saparam (void)
{
    GTNEXT();
    if (!strcmp(stream_->lexeme, "temp"))
        return COM_T;
    if (!strcmp(stream_->lexeme, "iter"))
        return COM_ITER;
    if (!strcmp(stream_->lexeme, "alpha"))
        return COM_ALPHA;
    if (!strcmp(stream_->lexeme, "beta"))
        return COM_BETA;
    if (!strcmp(stream_->lexeme, "perturb"))
        return COM_PERTURB;
    printf ("Expected 'temp', 'iter', 'alpha', or 'beta', but got %s\n", stream_->lexeme);
    return COM_ERR;
}

int sashow (void)
{
    int val;
    
    if (!strcmp(stream_->lexeme, "best")) {
        GTNEXT();
        val = p_feasible();
        if (val != COM_ERR)
            return 1 + val;
        return COM_ERR;
    }
    else if (stream_->type != T_EOF) {
        printf ("Expected 'best' or nothing, but got '%s'.\n", stream_->lexeme);
        return -1;
    }
    return 0;
}