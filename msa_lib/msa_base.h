#ifndef _MSA_BASE_H
#define _MSA_BASE_H

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <parasail.h>
#include <cassert>

// debug options
//#define DEBUG_MSA
//#define DEBUG_SEG
//#define DEBUG_MSA_FUNC

namespace msa_lib
{
#ifdef DEBUG_MSA
static const size_t LINE = 75;          // deprecated! line length of output
#endif
//static const size_t NUMBER = 200;       // deprecated! maximum number of sequences
static const size_t SIGMA = 128;        // alphabet size
static const char   DASH = '-';         // null symbol
#ifdef DEBUG_MSA
static const size_t BIG = 0x80000000;
#else
static const size_t BIG = 0x80000000;
#endif
static const unsigned DIAG = 0;         // code for traceback
static const unsigned VERT = 1;         // code for traceback
static const unsigned HORZ = 2;         // code for traceback
static const long long ROOT = -9;
static const long long INOD = -1;
static const long long MINE = 5;
static const long long MAXE = 50;


typedef struct coordinate COORDINATE;
typedef struct coordinate_values COORDINATE_VALUES;
typedef struct vertex VERTEX;
typedef struct edge EDGE;

// array for accessing vertices by lattice coordinate
struct coordinate{
    long long lo, hi;                         // lower and upper limits on array indices
    COORDINATE_VALUES *coord_vals;      // next coordinate array indices
    COORDINATE_VALUES *prev_coord_val;  // previous coordinate array index
    COORDINATE *next_on_free_list;      // maintains available records
    long long refer;                          // how many valid coordinate values do I have
};


// index in array
struct coordinate_values {
  COORDINATE	*next_coord;	// next coordinate array
  COORDINATE	*curr_coord;	// current coordinate array
  //long long value;          // value of this coordinate in absolute terms
};

// lattice vertex
struct vertex {
  EDGE *out;		                    // outgoing edge adjacency list
  COORDINATE_VALUES	*prev_coord_val;    // father in array
  EDGE *nonextracted_inedges;               // incoming edges still not extracted from heap
};

// lattice edge
struct edge {
	VERTEX	*tail, *head;	        // edge tail and head vertices
	long long	dist;		        // distance to head from source along edge
        long long     refer;                  // how many backtrack edges point to me
	EDGE *next, *prev;	        // edge adjacency list links
	EDGE *heap_succ, *heap_pred;	// heap bucket links
	EDGE *nonextracted_next,
             *nonextracted_prev;        // nonextracted_inedges links
	EDGE *backtrack;	        // edge to previous edge in path
};

// discrete heap of edges ranked by distance
class HEAP 
{
public:
	size_t min, max;	// minimum and maximum buckets of heap
	std::vector<EDGE *> bucket;	// buckets of edges
};

/* row in region on face of lattice */
class ROW
{
public:
    std::vector<long long> column;
    long long width;
    ROW(): width(0)
    {}
};

/*****************************************************************
*** bias.c                                                     ***
*****************************************************************/
class NODE{
public:
    long long sqn;
    double ltop;
    double w;
    double W;
    double v;
    double V;
    class NODE *lt;
    class NODE *rt;
    class NODE *par;
    class NODE *bro;

    NODE(): sqn(0), ltop(0), w(0), W(0), v(0), V(0), 
            lt(NULL), rt(NULL), par(NULL), bro(NULL)
    {}
};

/*****************************************************************
*** NDArray                                                    ***
*****************************************************************/
template<typename T>
class NDArray
{
public:
    int dim;
    std::vector<long long> dimens;
    std::vector<T> a;
public:
    NDArray(int dimension): dim(dimension)
    {}

    NDArray(int dimension, const std::vector<long long>& dimensions):
        dim(dimension), dimens( dimensions )
    {}

    void set_dimensions(long long m, long long n)
    {
        dimens.clear();
        dimens.push_back(m);
        dimens.push_back(n);
	a.assign(m * n, T());
    }

    void set_dimensions(long long m, long long n, long long p)
    {
        dimens.clear();
        dimens.push_back(m);
        dimens.push_back(n);
        dimens.push_back(p);
	a.assign(m * n * p, T());
    }

    T& at(long long i, long long j)
    {
    #ifdef DEBUG_SEG
        return a.at( i * dimens[1] + j );
    #else
        return a[ i * dimens[1] + j ];
    #endif
    }

    T& at(long long i, long long j, long long k)
    {
    #ifdef DEBUG_SEG
        return a.at( i * dimens[1] * dimens[2] + j * dimens[2] + k );
    #else
        return a[ i * dimens[1] * dimens[2] + j * dimens[2] + k ];
    #endif
    }
};

// class declaration
class MSA;
class MSA_Bias;
class MSA_Ecalc;
}// end of namespace msa_lib

#endif
