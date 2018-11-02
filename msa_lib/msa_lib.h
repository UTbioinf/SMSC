#ifndef _MSA_LIB_H
#define _MSA_LIB_H

#include "msa_base.h"
#include "msa_bias.h"
#include "msa_ecalc.h"

namespace msa_lib
{
/*****************************************************************
*** MSA                                                        ***
*****************************************************************/
class MSA
{
private:
    friend class MSA_Bias;
    friend class MSA_Ecalc;
    /*****************************************************************
    *** main.c                                                     ***
    *****************************************************************/

    //size_t NUMBER;    // NUMBER == K
    std::vector<std::string>  S;  // sequences: [NUMBER+1]
    long long delta;          // upper and lower bound differences
    NDArray<long long> epsi;  // projected & pairwise cost diff:  [NUMBER+1][NUMBER+1]
    NDArray<long long> scale; // pairwise cost weight scale: [NUMBER+1][NUMBER+1]
    NDArray<long long> Con;   // consistency check: [NUMBER+1][NUMBER+1][LENGTH+1]
    NDArray<long long> proj;  // projectd costs: [NUMBER+1][NUMBER+1]
    long long Upper, Lower;   // Upper and lower bounds on alignment distance
    NDArray<long long> dd,    // forward diagonal distance   [LENGTH+1]
                 hh,    // forward horizontal distance [LENGTH+1]
                 vv;    // forward vertical distance   [LENGTH+1]
    VERTEX *presource;  // Vertex before source; tail of first edge

    short D[SIGMA][SIGMA];      // symbol distance
    long long T[3][3][3][3];          // Altschul gap counts
    NDArray<long long *> Tpointer;    // Cuts off first two dimensions of T
    long long G, GG;                  // gap cost
    size_t max_len;
    HEAP h;

    NDArray<ROW> face;      // faces of lattice [NUMBER+1][NUMBER+1][LENGTH+1]
    NDArray<long long> costs;     // pairwise costs   [NUMBER+1][NUMBER+1]
    std::vector<COORDINATE *> msa_A;    // array to access vertices by lattice coordinate

    std::vector<std::string> M; // Multiple sequence alignment matrix
    long long C;                      // Current column in alignment matrix
    VERTEX *avail_vertex;
    EDGE *avail_edge;
    COORDINATE *avail_coordinate;

    std::vector<VERTEX* > vertex_pointers;
    std::vector<EDGE* > edge_pointers;
    std::vector<COORDINATE *> coordinate_pointers;

    EDGE* msa_result;
public:
    MSA();
    bool is_match_or_mismatch(size_t loc_s1, size_t loc_s2, const int s2Len,
            int gap_open, int gap_extend, 
            int match_score, int cur_score, const int* score_table);
    bool is_insertion(size_t loc_s1, size_t loc_s2, const int s2Len,
            int gap_open, int gap_extend,
            int cur_score, const int* score_table);
    bool is_deletion( size_t loc_s1, size_t loc_s2, const int s2Len,
            int gap_open, int gap_extend,
            int cur_score, const int* score_table);
    void nw_alignment(const std::string& s1, const std::string& s2, 
            int match, int mismatch, int gap_open, int gap_extend);
    void init();
    void set_NUC_matrix(long long match, long long mismatch, long long gap);
    void add_seq(const std::string& str);
    void run_msa(bool optimal);
    void free_msa();
    std::vector<std::string>& get_alignment();
    /*****************************************************************
    *** main.c                                                     ***
    *****************************************************************/
    long long min3(long long a, long long b, long long c);
    void INSERT(EDGE* &e, HEAP &h);
    void DELETE(EDGE* &e, HEAP &h);
    EDGE* msa();
    VERTEX* source();
    VERTEX* sink();
    long long intersect(std::vector<long long>& point, long long seqnum, std::vector<long long>& possible_values);
    void adjacent(EDGE* e, std::vector<long long>& q);
    void safe_coord(VERTEX* v, std::vector<long long>& p);
    void coord(VERTEX* v, std::vector<long long>& p);
    void project(std::vector<long long>& p, std::vector<long long>& q, std::vector<long long>& r);
    void compute_opt_aln(EDGE* e);
#ifdef DEBUG_MSA
    void display(EDGE* e);
    void column(std::vector<long long>& p, std::vector<long long>& q);
    void output();
#endif
    void heap(long long max);
    EDGE* extract();
    VERTEX* create_vertex(COORDINATE_VALUES* prev_coord_val);
    void free_vertex(VERTEX* v);
    EDGE* create_edge(VERTEX *v, VERTEX *w);
    void free_edge(EDGE* e);
    COORDINATE *create_coordinate(std::vector<long long>& index, long long n, COORDINATE_VALUES *prev_coord_val);
    bool free_coordinate(COORDINATE *a);

    void bias();
    void ecalc(long long len, bool optimal);
    /*****************************************************************
    *** primer.c                                                   ***
    *****************************************************************/
    void primer(); // may be replaced by Parasail later
    long long convert(long long I, long long J, long long n, long long m); // Calculate distance measure for evolutionary tree
    /*****************************************************************
    *** faces.c                                                    ***
    *****************************************************************/
    void faces();
};

}// namespace msa_lib end

#endif
