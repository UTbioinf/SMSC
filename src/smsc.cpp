#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <list>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <errno.h>

#include "common-bio/Progress.h"
#include "common-bio/Consensus.h"
#include "common-bio/Aligner.h"
#include "common-bio/Multiseq.h"
#include "common-bio/util.h"
#include "common-bio/baseheader.h"

#include "nucmer.h"

#include "mummer_api.h"

using namespace std;

#ifdef DEBUG_CHGCNT
size_t chgcnt = 0;
#endif

/* constants start */
/*coverage is high*/
const int MIN_OVERLAP = 50;
const double min_score = 0.99;
const double alpha = 0.9;

/*coverage is low*/
// const int MIN_OVERLAP = 1;
// const double min_score = 0.8;
// const double alpha = 0.7;

const double beta = (1-alpha) * 2;
/* constants end */

string ref_file, qry_file;
string filter_prefix = "";
string prefix = "smsc";
string ctg_file, idx_file;
int run_options = 4;
int consensus_strategy = 1;
int n_cores = 1;
ofstream fout_log;

loon::Progress progress;

size_t str2lu(const string& str)
{
    size_t ret = 0;
    for(size_t i = 0; i < str.length(); ++i)
        ret = ret * 10 + str[i] - '0';
    return ret;
}

class SegInfo
{
public:
    size_t ref_id, start_loc, end_loc;
public:
    SegInfo(size_t r_id = 0, size_t s_loc = 1, size_t e_loc = 0):
        ref_id(r_id), start_loc(s_loc), end_loc(e_loc)
    {
        if(s_loc == 0)
        {
            ERROR_PRINT("Implementation Error: start_loc starts from 0, not from 1");
            exit(1);
        }
    }

    void set_ref_id(size_t r_id)
    {
        ref_id = r_id;
    }

    void set_start_loc(size_t s_loc)
    {
        start_loc = s_loc;
    }

    void set_end_loc(size_t e_loc)
    {
        end_loc = e_loc;
    }
};

vector<SegInfo> seginfo_vec;

class Consensus: public loon::Consensus
{
public:
    list<size_t> path_forward;
    list<size_t> path_backward;
public:
    Consensus()
    {
    }

    Consensus(const loon::Consensus& base): loon::Consensus(base)
    {
    }

    void reverse_path()
    {
        path_forward.swap( path_backward );
    }

    void append_segment(size_t seg_id)
    {
        seg_id <<= 1;
        path_forward.push_back( seg_id );
        path_backward.push_back( seg_id + 1 );
    }

    void append_forward(Consensus& cons)
    {
        path_forward.splice(path_forward.end(), cons.path_forward);
        path_backward.splice(path_backward.begin(), cons.path_backward);
    }

    void append_backward(Consensus& cons)
    {
        cons.reverse_path();
        path_forward.splice(path_forward.end(), cons.path_forward);
        path_backward.splice(path_backward.begin(), cons.path_backward);
    }

    void swap(Consensus& cons)
    {
        path_forward.swap( cons.path_forward );
        path_backward.swap( cons.path_backward );
        loon::Consensus::swap( cons );
    }

    void move_forward_path(Consensus& cons)
    {
        path_forward.swap( cons.path_forward );
        path_backward.swap( cons.path_backward );
    }

    void move_backward_path(Consensus& cons)
    {
        path_forward.swap( cons.path_backward );
        path_backward.swap( cons.path_forward );
    }

    void print_path(ostream& out, size_t id)
    {
        out << id << ' ' << path_forward.size() << endl;
        for(list<size_t>::iterator p_it = path_forward.begin(); p_it != path_forward.end(); ++p_it)
        {
            size_t seg_id = *p_it;
            bool is_backward = seg_id & 1;
            seg_id >>= 1;
            out << seginfo_vec[ seg_id ].ref_id << ' '
                << seginfo_vec[ seg_id ].start_loc << ' '
                << seginfo_vec[ seg_id ].end_loc << ' '
                << is_backward << endl;
        }
    }
};

class KM
{
public:
    long long n;
    long long inf_val;
    long long res;
    loon::Array2D<long long> w;
    vector<bool> vis_x, vis_y;
    vector<long long> lx, ly;
    vector<long long> link_y;
    vector<long long> slack;
public:
    KM(): n(0), inf_val(0x7fffffff)
    {}
    void resize();
    void run();
    bool find(long long x);
};

void KM::resize()
{
    w.resize(n, n);
}

void KM::run()
{
    if(n == 0)
        return;

    long long i, j;
    link_y.assign(n, -1); // resize if necessary
    ly.assign(n, 0); // resize if necessary
    lx.assign(n, 0); // resize if necessary

    for(i = 0; i < n; ++i)
    {
        lx[i] = w.at(i, 0);
        for(j = 1; j < n; ++j)
            if(w.at(i, j) > lx[i])
                lx[i] = w.at(i, j);
        inf_val = max(inf_val, lx[i]);
    }

    for(long long x = 0; x < n; ++x)
    {
        slack.assign(n, inf_val); // resize if necessary

        while(true)
        {
            vis_x.assign(n, false); // resize if necessary
            vis_y.assign(n, false); // resize if necessary
            if(find(x))
                break;
            long long d = inf_val;
            for(i = 0; i < n; ++i)
                if(!vis_y[i] && d > slack[i])
                    d = slack[i];
            for(i = 0; i < n; ++i)
            {
                if(vis_x[i])
                    lx[i] -= d;

                if(vis_y[i])
                    ly[i] += d;
                else
                    slack[i] -= d;
            }
        }
    }

    res = 0;
    for(i = 0; i < n; ++i)
        if(link_y[i] > -1)
            res += w.at(link_y[ i ], i );
}

bool KM::find(long long x)
{
    vis_x[x] = true;
    for(int y = 0; y < n; ++y)
    {
        if(vis_y[y])
            continue;
        long long t = lx[x] + ly[y] - w.at(x, y);
        if(t == 0)
        {
            vis_y[y] = true;
            if(link_y[y] == -1 || find(link_y[y]))
            {
                link_y[y] = x;
                return true;
            }
        }
        else if(slack[y] > t)
            slack[y] = t;
    }
    return false;
}

/*======================= class ComEdge and ComEdges ============================
================================== start ======================================*/
class ComEdge
{
public:
#ifdef DEBUG
    static long long cnt;
#endif
    loon::Consensus comedge;
    set<size_t> left_vertices, right_vertices;
    long long weight;
    const loon::Fasta* fasta_stamp;
    bool stamp_is_forward;
public:
    ComEdge(): weight(0), fasta_stamp(NULL), stamp_is_forward(false)
    {}
    void consensus_edge(const loon::Fasta& fasta, const string& str, bool is_forward,
            vector<size_t>& left_vertices, vector<size_t>& right_vertices,
            vector<size_t>& left_difference, vector<size_t>& right_difference);
};



#ifdef DEBUG
long long ComEdge::cnt = 0;
#endif

void ComEdge::consensus_edge(const loon::Fasta& fasta, const string& str, bool is_forward,
        vector<size_t>& left_vertices, vector<size_t>& right_vertices,
        vector<size_t>& left_difference, vector<size_t>& right_difference)
{
    if(fasta_stamp == &fasta && stamp_is_forward == is_forward)
        throw 0x100;

    const double threshold = 0.85;
    vector<size_t>  t_vertices(this->left_vertices.begin(), this->left_vertices.end());
    set_difference(t_vertices.begin(), t_vertices.end(), left_vertices.begin(), left_vertices.end(), 
            inserter(left_difference, left_difference.begin()));
    size_t int_cnt = t_vertices.size() - left_difference.size();
    double union_cnt = t_vertices.size() + left_vertices.size() - int_cnt;
    if(int_cnt / union_cnt < threshold)
        throw 0x200;

    t_vertices.assign(this->right_vertices.begin(), this->right_vertices.end());
    set_difference(t_vertices.begin(), t_vertices.end(), right_vertices.begin(), right_vertices.end(),
            inserter(right_difference, right_difference.begin()));
    int_cnt = t_vertices.size() - right_difference.size();
    union_cnt = t_vertices.size() + right_vertices.size() - int_cnt;
    if(int_cnt / union_cnt < threshold)
        throw 0x300;

    fasta_stamp = &fasta;
    stamp_is_forward = is_forward;

#ifdef DEBUG
    ++cnt;
#endif
    if(comedge.size() == 0)
    {
        ERROR_PRINT("bazinga");
    }
    if(run_options & 0x10)
        comedge.consensus_full(loon::blasr_align_choice, str);
    else
        comedge.consensus_full( str );
    
    ++weight;
    for(vector<size_t>::iterator it = left_difference.begin(); it != left_difference.end(); ++it)
        this->left_vertices.erase( *it );
    for(vector<size_t>::iterator it = right_difference.begin(); it != right_difference.end(); ++it)
        this->right_vertices.erase( *it );
}

//------------------------ This is the dividing line ----------------------------

class SimpleVertex
{
public:
    set<size_t> succ_edges, pred_edges;
};

//------------------------ This is the dividing line ----------------------------

class PrimalGraph
{
public:
    vector<ComEdge> com_edges;
    vector<SimpleVertex> simp_vertices;
public:
    PrimalGraph(size_t n = 0): simp_vertices(n)
    {}
    void add_edge(const loon::Fasta& fasta, const string& str, bool is_forward, 
            vector<size_t>& left_vertices, vector<size_t>& right_vertices);
    void add_succ_edges(vector<size_t>& left_vertices, size_t e_id);
    void add_pred_edges(vector<size_t>& right_vertices, size_t e_id);
    void remove_succ_edges(vector<size_t>& left_rm, size_t e_id);
    void remove_pred_edges(vector<size_t>& right_rm, size_t e_id);
};

void PrimalGraph::add_edge(const loon::Fasta& fasta, const string& str, bool is_forward,
        vector<size_t>& left_vertices, vector<size_t>& right_vertices)
{
    for(size_t u_p = 0; u_p < left_vertices.size(); ++u_p)
    {
        size_t u_id = left_vertices[ u_p ];
        for(set<size_t>::iterator ep_it = simp_vertices[ u_id ].succ_edges.begin();
                ep_it != simp_vertices[ u_id ].succ_edges.end(); ++ep_it)
        {
            size_t e_id = *ep_it;
            try
            {
                vector<size_t> left_diff, right_diff;
                com_edges[e_id].consensus_edge(fasta, str, is_forward,
                        left_vertices, right_vertices,
                        left_diff, right_diff);
                remove_succ_edges(left_diff, e_id);
                remove_pred_edges(right_diff, e_id);
                return;
            }
            catch(int e)
            {
                continue;
            }
        }
    }
    for(size_t v_p = 0; v_p < right_vertices.size(); ++v_p)
    {
        size_t v_id = right_vertices[ v_p ];
        for(set<size_t>::iterator ep_it = simp_vertices[ v_id ].pred_edges.begin();
                ep_it != simp_vertices[ v_id ].pred_edges.end(); ++ep_it)
        {
            size_t e_id = *ep_it;
            try
            {
                vector<size_t> left_diff, right_diff;
                com_edges[e_id].consensus_edge(fasta, str, is_forward,
                        left_vertices, right_vertices,
                        left_diff, right_diff);
                remove_succ_edges(left_diff, e_id);
                remove_pred_edges(right_diff, e_id);
            }
            catch(int e)
            {
                continue;
            }
        }
    }

    // new edge
    size_t e_id = com_edges.size();
    com_edges.push_back( ComEdge() );
    com_edges[ e_id ].comedge.set_str( str );
    com_edges[ e_id ].weight = 1;
    com_edges[ e_id ].fasta_stamp = &fasta;
    com_edges[ e_id ].stamp_is_forward = is_forward;
    com_edges[ e_id ].left_vertices.insert(left_vertices.begin(), left_vertices.end() );
    com_edges[ e_id ].right_vertices.insert(right_vertices.begin(), right_vertices.end() );

    add_succ_edges(left_vertices, e_id);
    add_pred_edges(right_vertices, e_id);
}

void PrimalGraph::add_succ_edges(vector<size_t>& left_vertices, size_t e_id)
{
    for(size_t u_p = 0; u_p < left_vertices.size(); ++u_p)
        simp_vertices[ left_vertices[u_p] ].succ_edges.insert( e_id );
}

void PrimalGraph::add_pred_edges(vector<size_t>& right_vertices, size_t e_id)
{
    for(size_t v_p = 0; v_p < right_vertices.size(); ++v_p)
        simp_vertices[ right_vertices[ v_p ] ].pred_edges.insert( e_id );
}

void PrimalGraph::remove_succ_edges(vector<size_t>& left_rm, size_t e_id)
{
    for(vector<size_t>::iterator it = left_rm.begin(); it != left_rm.end(); ++it)
        simp_vertices[ *it ].succ_edges.erase( e_id );
}

void PrimalGraph::remove_pred_edges(vector<size_t>& right_rm, size_t e_id)
{
    for(vector<size_t>::iterator it = right_rm.begin(); it != right_rm.end(); ++it)
        simp_vertices[ *it ].pred_edges.erase( e_id );
}
/*====================== class ComEdge and ComEdges ============================
================================== end =======================================*/

/*============================== class CodirGraph ==============================
====================================== start =================================*/
class CodirGraph
{
public:
    vector<size_t> vertex_xtofx;
    vector<size_t> unused_vertices;
    vector<size_t> vertex_fxtox;
    vector<map<size_t, ComEdge*> > graph;
    KM km;
    vector<vector<size_t> > paths;
public:
    CodirGraph(PrimalGraph& p_graph, vector<Consensus>& contigs);
    void run();
    void get_path_head(size_t i, size_t& u);
    void get_edge_and_vertex(size_t i, size_t j, ComEdge* &e, size_t& v);
};

CodirGraph::CodirGraph(PrimalGraph& p_graph, vector<Consensus>& contigs)
{
    vertex_xtofx.resize( p_graph.simp_vertices.size(), loon::INF );
    for(size_t x = 0; x < p_graph.simp_vertices.size(); ++x)
    {
        if(!p_graph.simp_vertices[x].succ_edges.empty() || 
            !p_graph.simp_vertices[x].pred_edges.empty())
        {
            vertex_xtofx[x] = vertex_fxtox.size();
            vertex_fxtox.push_back( x );
        }
        else
        {
            unused_vertices.push_back( x );
        }
    }

    size_t n = vertex_fxtox.size();
    graph.resize( n );

    loon::Aligner verify_aligner;
    verify_aligner.set_simple_align();
    for(size_t fx = 0; fx < n; ++fx)
    {
        size_t x = vertex_fxtox[ fx ];
        for(set<size_t>::iterator succ_e_it = p_graph.simp_vertices[ x ].succ_edges.begin();
                succ_e_it != p_graph.simp_vertices[ x ].succ_edges.end(); ++succ_e_it)
        {
            ComEdge* e = &p_graph.com_edges[ *succ_e_it ];
            // verify before considering this edge
            verify_aligner.set_strs( contigs[ x ].toString(), e->comedge.toString() );
            if(!((run_options & 0x10) ? 
                    verify_aligner.align_right(loon::blasr_align_choice, true) : 
                    verify_aligner.align_right(true)))
                continue;

            for(set<size_t>::iterator vx_it = e->right_vertices.begin();
                    vx_it != e->right_vertices.end(); ++vx_it)
            {
                if(x == *vx_it)
                    continue;
                size_t vfx = vertex_xtofx[ *vx_it ];
                map<size_t, ComEdge*>::iterator succ_it = graph[ fx ].find( vfx );
                if(succ_it == graph[fx].end())
                {
                    // verify before adding the edge
                    verify_aligner.set_strs( contigs[ *vx_it ].toString(), e->comedge.toString() );
                    if(!((run_options & 0x10) ?
                            verify_aligner.align_left(loon::blasr_align_choice, true) : 
                            verify_aligner.align_left(true)))
                        continue;
                    graph[fx][ vfx ] = e;
                }
                else
                {
                    if(succ_it->second->weight < e->weight)
                    {
                        // verify before adding the edge
                        verify_aligner.set_strs( contigs[  *vx_it ].toString(), e->comedge.toString() );
                        if(!((run_options & 0x10) ?
                                verify_aligner.align_left(loon::blasr_align_choice, true) : 
                                verify_aligner.align_left(true)))
                            continue;
                        graph[fx][vfx] = e;
                    }
                }
            }
        }
    }

    km.n = n;
    km.resize();
    km.w.assign(n * n, 0);

    for(size_t ufx = 0; ufx < n; ++ufx)
        for(map<size_t, ComEdge*>::iterator succ_it = graph[ufx].begin();
                succ_it != graph[ ufx ].end(); ++succ_it)
            km.w.at(ufx, succ_it->first) = succ_it->second->weight;
}

void CodirGraph::run()
{
    km.run();

    vector<size_t> link_x( km.n );
    vector<bool> visited_x( km.n, false );
    for(size_t i = 0; i < km.n; ++i)
        link_x[ km.link_y[i] ] = i;

    for(size_t i = 0; i < link_x.size(); ++i)
    {
        if(!visited_x[i])
        {
            vector<size_t> fp, bp;
            // find forward
            bool is_cycle = false;
            size_t u = i;
            size_t min_edge_tail;

            visited_x[ u ] = true;
            while( km.w.at(u, link_x[u]) > 0)
            {
                if(fp.empty())
                    min_edge_tail = 0;
                else if(km.w.at( km.link_y[ fp[min_edge_tail] ] , fp[min_edge_tail] ) > km.w.at( u , link_x[u] ))
                    min_edge_tail = fp.size();

                u = link_x[ u ];
                fp.push_back( u );


                if(visited_x[ u ])
                {
                    is_cycle = true;
                    break;
                }

                visited_x[ u ] = true;
            }

            if(is_cycle)
            {
                paths.push_back( vector<size_t>() );
                vector<size_t>& one_path = paths.back();
                one_path.assign( fp.begin() + min_edge_tail, fp.end() );
                one_path.insert( one_path.end(),  fp.begin(), fp.begin() + min_edge_tail );
                continue;
            }

            // find backward
            u = i;
            while(km.w.at( km.link_y[u] , u) > 0)
            {
                u = km.link_y[u];
                bp.push_back( u );
                visited_x[ u ] = true;
            }

            paths.push_back( vector<size_t>() );
            vector<size_t>& one_path = paths.back();
            one_path.assign(bp.rbegin(), bp.rend());
            one_path.push_back(i);
            one_path.insert(one_path.end(), fp.begin(), fp.end());
        }
    }

    DEBUG_PRINT("output paths begin: n = %lld", km.n);
#ifdef DEBUG
    for(size_t i = 0; i < paths.size(); ++i)
    {
        cerr << "Path " << i << ": " << paths[i][0];
        for(size_t j = 1; j < paths[i].size(); ++j)
                cerr << " -> " << paths[i][j];
        cerr << endl;
    }
#endif
    DEBUG_PRINT("output paths end");
}

void CodirGraph::get_path_head(size_t i, size_t& u)
{
    u = vertex_fxtox[ paths[i][0] ];
}

void CodirGraph::get_edge_and_vertex(size_t i, size_t j, ComEdge* &e, size_t& v)
{
    e = graph[ paths[i][j-1] ][ paths[i][j] ];
    v = vertex_fxtox[ paths[i][j] ];
}
/*============================== class CodirGraph ==============================
======================================= end ==================================*/

bool add_edge(vector<size_t>& left_vertices, vector<size_t>& right_vertices,
        const tigrinc::Nucmer_Delta& delta)
{
    long long ref_left_ext = delta.clusters.front().ref_start - 1;
    long long ref_right_ext = delta.ref_len - delta.clusters.back().ref_end;
    long long qry_left_ext = delta.clusters.front().qry_start - 1;
    long long qry_right_ext = delta.qry_len - delta.clusters.back().qry_end;

    if(ref_left_ext >= qry_left_ext && ref_right_ext >= qry_right_ext)
        return false;
    if(ref_left_ext <= qry_left_ext && ref_right_ext <= qry_right_ext)
        return false;

    if(qry_left_ext > ref_left_ext)
        right_vertices.push_back( str2lu( delta.ref_tag ) );
    else
        left_vertices.push_back( str2lu(delta.ref_tag) );
    return true;
}

void bridge_codir(vector<Consensus>& contigs,
        loon::Multiseq& long_reads, 
        vector<Consensus>& codir_res,
        vector<bool>& used_long_reads)
{
    const size_t lr_n = long_reads.size();
    used_long_reads.assign(lr_n, false);

    PrimalGraph p_graph( contigs.size() );

    loon::Aligner aligner;
    aligner.set_simple_align();

    INFO_PRINT("Convert contigs to strings");
    loon::Multiseq tmp_ref;
    tmp_ref.resize( contigs.size() );
    for(size_t j = 0; j < contigs.size(); ++j)
        tmp_ref[j].seq = contigs[j].toString();

    INFO_PRINT("All contigs vs. all long reads matching");
    nucmer_Deltas deltas_forward, deltas_backward;
    
    loon::all_v_all_matching(tmp_ref, long_reads, deltas_forward, deltas_backward, true);

    vector<vector<size_t> > fleft_vertices( lr_n ), fright_vertices( lr_n );
    vector<vector<size_t> > rleft_vertices( lr_n ), rright_vertices( lr_n );

    // compute vertices for edges
    INFO_PRINT("Compute vertices for edges");
    progress.start( deltas_forward.size() + deltas_backward.size() );
    for(size_t ii = 0; ii < deltas_forward.size(); ++ii)
    {
        size_t lr_i = str2lu(deltas_forward[ii].qry_tag);
        used_long_reads[ lr_i ] = add_edge(fleft_vertices[ lr_i ], 
                                            fright_vertices[ lr_i ], 
                                            deltas_forward[ ii ]) || 
                                    used_long_reads[ lr_i ];

        progress.progress(1);
    }

    for(size_t ii = 0; ii < deltas_backward.size(); ++ii)
    {
        size_t lr_i = str2lu(deltas_backward[ ii ].qry_tag);
        used_long_reads[ lr_i ] = add_edge(rleft_vertices[ lr_i ], 
                                            rright_vertices[ lr_i ], 
                                            deltas_backward[ ii ]) || 
                                    used_long_reads[ lr_i ];

        progress.progress(1);
    }
    progress.stop();

    // build the primal graph
    INFO_PRINT("Build the primal graph");

    progress.start( lr_n );
    for(size_t lr_i = 0; lr_i < lr_n; ++lr_i)
    {
        if(!fleft_vertices[ lr_i ].empty() && !fright_vertices[ lr_i ].empty())
        {
            p_graph.add_edge( long_reads[ lr_i ], long_reads[ lr_i ].seq, true, 
                    fleft_vertices[ lr_i ], fright_vertices[ lr_i ]);
        }
        if(!rleft_vertices[ lr_i ].empty() && !rright_vertices[ lr_i ].empty())
        {
            p_graph.add_edge( long_reads[ lr_i ], long_reads[ lr_i ].rev, true,
                    rleft_vertices[ lr_i ], rright_vertices[ lr_i ]);
        }
        progress.progress(1);
    }

    DEBUG_PRINT("Print primal graph");
    DEBUG_PRINT(">>> edge oriented: %lu edges", p_graph.com_edges.size());
#ifdef DEBUG
    for(size_t i = 0; i < p_graph.com_edges.size(); ++i)
    {
        DEBUG_PRINT("edge id = %lu, weight = %lld", i, p_graph.com_edges[i].weight);
        DEBUG_PRINT("left vertices:");
        for(set<size_t>::iterator lv_it = p_graph.com_edges[i].left_vertices.begin();
                lv_it != p_graph.com_edges[i].left_vertices.end(); ++lv_it)
            cerr << (*lv_it) << ' ';
        cerr << endl;
        DEBUG_PRINT("right vertices:");
        for(set<size_t>::iterator rv_it = p_graph.com_edges[i].right_vertices.begin();
                rv_it != p_graph.com_edges[i].right_vertices.end(); ++rv_it)
            cerr << (*rv_it) << ' ';
        cerr << endl;
    }

    DEBUG_PRINT(">>> vertex oriented: %lu vertices", p_graph.simp_vertices.size());
    for(size_t i = 0; i < p_graph.simp_vertices.size(); ++i)
    {
        DEBUG_PRINT("vertex id = %lu", i);
        DEBUG_PRINT("pred edges (edge id, edge weight):");
        for(set<size_t>::iterator pe_it = p_graph.simp_vertices[i].pred_edges.begin();
                pe_it != p_graph.simp_vertices[i].pred_edges.end(); ++pe_it)
            cerr << "(" << (*pe_it) << ", " << p_graph.com_edges[ *pe_it ].weight << "); ";
        cerr << endl;

        DEBUG_PRINT("succ edges (edge id, edge weight):");
        for(set<size_t>::iterator se_it = p_graph.simp_vertices[i].succ_edges.begin();
                se_it != p_graph.simp_vertices[i].succ_edges.end(); ++se_it)
            cerr << "(" << (*se_it) << ", " << p_graph.com_edges[ *se_it ].weight << ") ";
        cerr << endl;
    }
#endif

    INFO_PRINT("Build the actual graph based on the primal graph");
    // build the actual graph based on the primal graph
    CodirGraph codir_graph(p_graph, contigs);

    DEBUG_PRINT("Print actual graph");
    DEBUG_PRINT("Vertex translation (fx to x)");
#ifdef DEBUG
    for(size_t fx = 0; fx < codir_graph.vertex_fxtox.size(); ++fx)
        DEBUG_PRINT("fx = %lu\tx = %lu", fx, codir_graph.vertex_fxtox[ fx ]);
    DEBUG_PRINT("The graph");
    for(size_t i = 0; i < codir_graph.graph.size(); ++i)
    {
        DEBUG_PRINT("vertex id: %lu", i);
        for(map<size_t, ComEdge*>::iterator v_it = codir_graph.graph[i].begin();
                v_it != codir_graph.graph[i].end(); ++v_it)
        {
            cerr << "(" << (v_it->first) << ", " << (v_it->second->weight) << "); ";
        }
        cerr << endl;
    }
#endif

    INFO_PRINT("Run max path cover (KM) algorithm");
    codir_graph.run();

    // do consensus
    INFO_PRINT("Do consensus");
    size_t u_id, v_id;
    ComEdge* e;
    for(size_t i = 0; i < codir_graph.paths.size(); ++i)
    {
        codir_graph.get_path_head(i, u_id);
        codir_res.push_back( Consensus() );
        codir_res.back().swap( contigs[u_id] );
        
        for(size_t j = 1; j < codir_graph.paths[i].size(); ++j)
        {
            codir_graph.get_edge_and_vertex(i, j, e, v_id);
            size_t edge_length = e->comedge.size();

            try
            {
                DEBUG_PRINT("concatenate left");
                if(run_options & 0x10)
                    codir_res.back().consensus_right( loon::blasr_align_choice, e->comedge.toString() );
                else
                    codir_res.back().consensus_right( e->comedge.toString() );
            }
            catch(int error_code)
            {
                DEBUG_PRINT("consensus_right error [1]: error_code = %d", error_code);
                codir_res.push_back( Consensus(e->comedge) );
            }

            try
            {
                DEBUG_PRINT("concatenate right");
                if(run_options & 0x10)
                {
                    codir_res.back().consensus_right(loon::blasr_align_choice, 
                                                    contigs[v_id].toString(), 
                                                    codir_res.back().size() - min(edge_length, 
                                                    contigs[v_id].size()) );
                }
                else
                {
                    codir_res.back().consensus_right( contigs[v_id].toString(), 
                                                        codir_res.back().size() - min(edge_length, 
                                                        contigs[v_id].size()) );
                }
                codir_res.back().append_forward( contigs[ v_id ] );

            }
            catch(int error_code)
            {
                DEBUG_PRINT("consensus_right error [2]: error_code = %d", error_code);
                codir_res.push_back( Consensus() );
                codir_res.back().swap( contigs[v_id] );
            }
        }
    }

    // collect unused vertices
    INFO_PRINT("Collect unused vertices");
    for(size_t i = 0; i < codir_graph.unused_vertices.size(); ++i)
    {
        codir_res.push_back( Consensus() );
        codir_res.back().swap( contigs[codir_graph.unused_vertices[i]] );
    }

    DEBUG_PRINT("ComEdge::cnt = %lld", ComEdge::cnt);
}

/*========================= class RPrimalGraph =============================
================================== start =================================*/
class RComEdge
{
public:
    loon::Consensus rcom_edge;
    set<size_t> type0_heads, type0_tails;
    set<size_t> type1_heads, type1_tails;
    long long weight;
    const loon::Fasta* fasta_stamp;
    bool stamp_is_forward;
public:
    RComEdge(): weight(0), fasta_stamp(NULL), stamp_is_forward(false)
    {}
    void consensus_edge(const loon::Fasta& fasta, const string& str, bool is_forward, bool is_type0,
            vector<size_t>& t_tails, vector<size_t>& t_heads,
            vector<size_t>& tails_diff, vector<size_t>& heads_diff);
    vector<size_t> set_vec_diff(const set<size_t>& s, const vector<size_t>& v);
};

void RComEdge::consensus_edge(const loon::Fasta& fasta, const string& str, bool is_forward, bool is_type0,
        vector<size_t>& t_tails, vector<size_t>& t_heads,
        vector<size_t>& tails_diff, vector<size_t>& heads_diff)
{
    if(fasta_stamp == &fasta && stamp_is_forward == is_forward)
        throw 0x100;

    const double threshold = 0.85;
    size_t int_cnt;
    double union_cnt;
    if(is_type0)
    {
        tails_diff = set_vec_diff( type0_tails, t_tails );
        int_cnt = type0_tails.size() - tails_diff.size();
        union_cnt = type0_tails.size() + t_tails.size() - int_cnt;
        if(int_cnt / union_cnt < threshold)
            throw 0x200;

        heads_diff = set_vec_diff( type0_heads, t_heads );
        int_cnt = type0_heads.size() - heads_diff.size();
        union_cnt = type0_heads.size() + t_heads.size() - int_cnt;
        if(int_cnt / union_cnt < threshold)
            throw 0x300;
    }
    else
    {
        tails_diff = set_vec_diff( type0_tails, t_tails );
        int_cnt = type1_tails.size() - tails_diff.size();
        union_cnt = type1_tails.size() + t_tails.size() - int_cnt;
        if(int_cnt / union_cnt < threshold)
            throw 0x200;

        heads_diff = set_vec_diff( type1_heads, t_heads );
        int_cnt = type1_heads.size() - heads_diff.size();
        union_cnt = type1_heads.size() + t_heads.size() - int_cnt;
        if(int_cnt / union_cnt < threshold)
            throw 0x300;
    }

    fasta_stamp = &fasta;
    stamp_is_forward = is_forward;

    if(run_options & 0x10)
        rcom_edge.consensus_full( loon::blasr_align_choice, str );
    else
        rcom_edge.consensus_full( str );
    ++weight;
}

vector<size_t> RComEdge::set_vec_diff(const set<size_t>& s, const vector<size_t>& v)
{
    vector<size_t> ret;
    set<size_t>::const_iterator s_it = s.begin();
    vector<size_t>::const_iterator v_it = v.begin();
    while(s_it != s.end() && v_it != v.end())
    {
        if(*s_it < *v_it)
        {
            ret.push_back( *s_it );
            ++s_it;
        }
        else if(*v_it < *s_it)
            ++v_it;
        else
        {
            ++s_it;
            ++v_it;
        }
    }

    ret.insert(ret.end(), s_it, s.end());
    return ret;
}

// ------------------------- This is the dividing line ---------------------------

class VertexPair
{
public:
    size_t u, v;
    bool is_type0;
public:
    VertexPair(size_t uu = loon::INF, size_t vv = loon::INF, bool type0 = false):
        u(uu), v(vv), is_type0(type0)
    {}

    friend bool operator<(const VertexPair& vp1, const VertexPair& vp2);
    friend ostream& operator<<(ostream& os, const VertexPair& vp);
};

bool operator<(const VertexPair& vp1, const VertexPair& vp2)
{
    return (vp1.u < vp2.u || (vp1.u == vp2.u && vp1.v < vp2.v) ||
            (vp1.u == vp2.u && vp1.v == vp2.v && vp1.is_type0 && !vp2.is_type0));
}

ostream& operator<<(ostream& os, const VertexPair& vp)
{
    os << "(" << vp.u << ", " << vp.v << ")";
    return os;
}

// ------------------------- This is the dividing line ---------------------------

class RPrimalGraph
{
public:
    vector<RComEdge> rcom_edges;
    map<VertexPair, set<size_t> > multiedges;
public:
    void add_type0_edge(const loon::Fasta& fasta, const string& str, const string& rev_str,
            vector<size_t>& fleft_vertices, vector<size_t>& rleft_vertices);
    void add_type1_edge(const loon::Fasta& fasta, const string& str, const string& rev_str,
            vector<size_t>& fright_vertices, vector<size_t>& rright_vertices);
    size_t get_new_edge_id(size_t edge_id, const string& str, 
            const loon::Fasta& fasta, 
            const bool is_forward = false,
            const long long weight = 1);
    void consensus_edge(set<size_t>& mes, const loon::Fasta& fasta,
            const string& str, const bool dir, size_t& edge_id,
            VertexPair& vp,
            vector<size_t>& t_tails, vector<size_t>& t_heads);
};

void RPrimalGraph::consensus_edge(set<size_t>& mes, const loon::Fasta& fasta, 
        const string& str, const bool dir, size_t& edge_id,
        VertexPair& vp,
        vector<size_t>& t_tails, vector<size_t>& t_heads)
{
    VertexPair erase_vp;
    erase_vp.is_type0 = vp.is_type0;

    for(set<size_t>::iterator meid_it = mes.begin(); meid_it != mes.end(); ++meid_it)
    {
        try
        {
            vector<size_t> tails_diff, heads_diff;
            rcom_edges[ *meid_it ].consensus_edge(fasta, str, dir, vp.is_type0,
                    t_tails, t_heads,
                    tails_diff, heads_diff);

            for(vector<size_t>::const_iterator t_it = tails_diff.begin(); t_it != tails_diff.end(); ++t_it)
                for(vector<size_t>::const_iterator h_it = heads_diff.begin(); h_it != heads_diff.end(); ++h_it)
                {
                    erase_vp.u = *t_it;
                    erase_vp.v = *h_it;
                    map<VertexPair, set<size_t> >::iterator map_it = multiedges.find( erase_vp );
                    if(map_it != multiedges.end())
                    {
                        map_it->second.erase( *meid_it );
                        if(map_it->second.empty())
                            multiedges.erase( map_it );
                    }
                }
            
            return;
        }
        catch(int error_code)
        {
            // continue;
        }
    }

    edge_id = get_new_edge_id( edge_id, str, fasta, dir, 1);
    mes.insert( edge_id );

    if(vp.is_type0)
    {
        rcom_edges[ edge_id ].type0_tails.insert( vp.u );
        rcom_edges[ edge_id ].type0_heads.insert( vp.v );
    }
    else
    {
        rcom_edges[ edge_id ].type1_tails.insert( vp.u );
        rcom_edges[ edge_id ].type1_heads.insert( vp.v );
    }
}

size_t RPrimalGraph::get_new_edge_id(size_t edge_id, const string& str,
        const loon::Fasta& fasta,
        const bool is_forward/* = false */,
        const long long weight/* = 1*/)
{
    if(edge_id == loon::INF)
    {
        edge_id = rcom_edges.size();
        rcom_edges.push_back( RComEdge() );
        rcom_edges[ edge_id ].rcom_edge.set_str( str );
        rcom_edges[ edge_id ].fasta_stamp = &fasta;
        rcom_edges[ edge_id ].stamp_is_forward = is_forward;
        rcom_edges[ edge_id ].weight = weight;
    }
    return edge_id;
}

void RPrimalGraph::add_type0_edge(const loon::Fasta& fasta, const string& str, const string& rev_str,
        vector<size_t>& fleft_vertices, vector<size_t>& rleft_vertices)
{
    VertexPair vp;
    vp.is_type0 = true;

    size_t new_fedge_id = loon::INF;
    size_t new_redge_id = loon::INF;

    for(vector<size_t>::const_iterator flv_it = fleft_vertices.begin(); 
            flv_it != fleft_vertices.end(); ++flv_it)
    {
        for(vector<size_t>::const_iterator rlv_it = rleft_vertices.begin();
                rlv_it != rleft_vertices.end(); ++rlv_it)
        {
            if(*flv_it == *rlv_it)
                continue;
            vp.u = *flv_it;
            vp.v = *rlv_it;

            map<VertexPair, set<size_t> >::iterator me_it = multiedges.find( vp );
            if(me_it != multiedges.end())
            { // try to consensus edge
                consensus_edge(me_it->second, fasta, str, true, new_fedge_id, vp,
                        fleft_vertices, rleft_vertices);
                continue;
            }

            swap(vp.u, vp.v);
            me_it = multiedges.find( vp );
            if(me_it != multiedges.end())
            { // try to consensus edge
                consensus_edge(me_it->second, fasta, rev_str, false, new_redge_id, vp,
                        rleft_vertices, fleft_vertices);
                continue;
            }

            swap(vp.v, vp.u);
            new_fedge_id = get_new_edge_id( new_fedge_id, str, fasta, true, 1);
            multiedges[ vp ] = set<size_t>();
            multiedges[ vp ].insert( new_fedge_id );

            rcom_edges[ new_fedge_id ].type0_tails.insert( vp.u );
            rcom_edges[ new_fedge_id ].type0_heads.insert( vp.v );
        }
    }
}

void RPrimalGraph::add_type1_edge(const loon::Fasta& fasta, const string& str, const string& rev_str,
        vector<size_t>& fright_vertices, vector<size_t>& rright_vertices)
{
    VertexPair vp;
    vp.is_type0 = false;

    size_t new_fedge_id = loon::INF;
    size_t new_redge_id = loon::INF;

    for(vector<size_t>::const_iterator rrv_it = rright_vertices.begin();
            rrv_it != rright_vertices.end(); ++rrv_it)
    {
        for(vector<size_t>::const_iterator frv_it = fright_vertices.begin();
                frv_it != fright_vertices.end(); ++frv_it)
        {
            if(*rrv_it == *frv_it)
                continue;
            vp.u = *rrv_it;
            vp.v = *frv_it;

            map<VertexPair, set<size_t> >::iterator me_it = multiedges.find( vp );
            if(me_it != multiedges.end())
            { // try to consensus edge
                consensus_edge(me_it->second, fasta, str, true, new_fedge_id, vp,
                        rright_vertices, fright_vertices);
                continue;
            }

            swap(vp.u, vp.v);
            me_it = multiedges.find( vp );
            if(me_it != multiedges.end())
            { // try to consensus edge
                consensus_edge(me_it->second, fasta, rev_str, false, new_redge_id, vp,
                        fright_vertices, rright_vertices);
                continue;
            }

            swap(vp.u, vp.v);
            new_fedge_id = get_new_edge_id( new_fedge_id, str, fasta, true, 1);
            multiedges[ vp ] = set<size_t>();
            multiedges[ vp ].insert( new_fedge_id );

            rcom_edges[ new_fedge_id ].type1_tails.insert( vp.u );
            rcom_edges[ new_fedge_id ].type1_heads.insert( vp.v );
        }
    }
}

/*========================= class RPrimalGraph ============================
=================================== end =================================*/

/*========================= class RevdirGraph =============================
================================= start =================================*/
class RevdirGraph
{
public:
    struct GraphEdge
    {
        RComEdge* rce_p;
        bool is_type0;

        GraphEdge(RComEdge* rp = NULL, bool type0 = false): rce_p(rp), is_type0(type0)
        {}
    };

    size_t n;
    vector<size_t> vertex_xtofx;
    vector<size_t> vertex_fxtox;
    vector<map<size_t, GraphEdge> > succ_edges;
    vector<map<size_t, GraphEdge> > pred_edges;
    KM km;
    vector<vector<size_t> > paths;
public:
    size_t get_vertex_original_id(size_t vid);
    size_t get_vertex_id(size_t v);
    bool find_edge(RComEdge* &ep, bool& is_type0, size_t uid, size_t vid); // return true means the edge direction should not be changed; return false otherwise
    void assign_vertex_id(size_t v);
    RevdirGraph(size_t nn, RPrimalGraph& rp_graph, vector<Consensus>& contigs);
    void run_max_path_cover();
};

size_t RevdirGraph::get_vertex_original_id(size_t vid)
{
    return vertex_fxtox[ vid ];
}

size_t RevdirGraph::get_vertex_id(size_t v)
{
    return vertex_xtofx[ v ];
}

bool RevdirGraph::find_edge(RComEdge* &ep, bool& is_type0, size_t uid, size_t vid)
{
    map<size_t, GraphEdge>::iterator ret_it = succ_edges[uid].find( vid );
    if(ret_it != succ_edges[uid].end())
    {
        ep = ret_it->second.rce_p;
        is_type0 = ret_it->second.is_type0;
        return true;
    }

    ret_it = pred_edges[uid].find( vid );
    if(ret_it != pred_edges[uid].end())
    {
        ep = ret_it->second.rce_p;
        is_type0 = ret_it->second.is_type0;
        return false;
    }
    
    throw 0x1;
}

void RevdirGraph::assign_vertex_id(size_t v)
{
    if(vertex_xtofx[ v ] == loon::INF)
    {
        vertex_xtofx[ v ] = vertex_fxtox.size();
        vertex_fxtox.push_back( v );

        succ_edges.push_back( map<size_t, GraphEdge>() );
        pred_edges.push_back( map<size_t, GraphEdge>() );
    }
}

RevdirGraph::RevdirGraph(size_t nn, RPrimalGraph& rp_graph, vector<Consensus>& contigs): n(nn)
{
    vertex_xtofx.resize(nn, loon::INF);

    GraphEdge tmp_edge;
    
    loon::Aligner verify_aligner;
    verify_aligner.set_simple_align();

    for(map<VertexPair, set<size_t> >::iterator me_it = rp_graph.multiedges.begin();
            me_it != rp_graph.multiedges.end(); ++me_it)
    {
        const VertexPair& vp = me_it->first;
        tmp_edge.rce_p = NULL;
        tmp_edge.is_type0 = vp.is_type0;
        for(set<size_t>::iterator meid_it = me_it->second.begin();
                meid_it != me_it->second.end(); ++meid_it)
        {
            if(tmp_edge.rce_p == NULL || tmp_edge.rce_p->weight < rp_graph.rcom_edges[ *meid_it ].weight)
            {
                string tt_str = rp_graph.rcom_edges[ *meid_it ].rcom_edge.toString();
                if(vp.is_type0)
                {
                    verify_aligner.set_strs(contigs[vp.u].toString(), tt_str);
                    if(!((run_options & 0x10) ? 
                            verify_aligner.align_right(loon::blasr_align_choice, true) : 
                            verify_aligner.align_right(true)))
                    {
                        DEBUG_PRINT_FILL_BLANK(1, "(%lu, %lu, 0): %lu -> (%lld) fails (error_code = %d)",
                                vp.u, vp.v, vp.u, rp_graph.rcom_edges[ *meid_it ].weight, verify_aligner.error_code);
                        continue;
                    }
                    loon::compute_reverse_complement(tt_str);
                    verify_aligner.set_strs(contigs[vp.v].toString(), tt_str);
                    if(!((run_options & 0x10) ?
                            verify_aligner.align_right(loon::blasr_align_choice, true) : 
                            verify_aligner.align_right(true)))
                    {
                        DEBUG_PRINT_FILL_BLANK(1, "(%lu, %lu, 0): (%lld) -> %lu\' fails (error_code = %d)",
                                vp.u, vp.v, rp_graph.rcom_edges[ *meid_it ].weight, vp.v, verify_aligner.error_code);
                        continue;
                    }
                }
                else
                {
                    verify_aligner.set_strs(contigs[vp.v].toString(), tt_str);
                    if(!((run_options & 0x10) ?
                            verify_aligner.align_left(loon::blasr_align_choice, true) :
                            verify_aligner.align_left(true) ))
                    {
                        DEBUG_PRINT_FILL_BLANK(1, "(%lu, %lu, 1): (%lld) -> %lu fails (error_code = %d)",
                                vp.u, vp.v, rp_graph.rcom_edges[ *meid_it ].weight, vp.v, verify_aligner.error_code);
                        continue;
                    }
                    loon::compute_reverse_complement(tt_str);
                    verify_aligner.set_strs(contigs[vp.u].toString(), tt_str);
                    if(!( (run_options & 0x10) ? 
                            verify_aligner.align_left(loon::blasr_align_choice, true) : 
                            verify_aligner.align_left(true)))
                    {
                        DEBUG_PRINT_FILL_BLANK(1, "(%lu, %lu, 1): %lu\' -> (%lld) fails (error_code = %d)",
                                vp.u, vp.v, vp.u, rp_graph.rcom_edges[ *meid_it ].weight, verify_aligner.error_code);
                        continue;
                    }
                }
                tmp_edge.rce_p = &rp_graph.rcom_edges[ *meid_it ];
            }
        }
        
        if(tmp_edge.rce_p == NULL)
            continue;

        assign_vertex_id( vp.u );
        assign_vertex_id( vp.v );

        size_t vpu_id = get_vertex_id(vp.u);
        size_t vpv_id = get_vertex_id(vp.v);
        map<size_t, GraphEdge>::iterator edge_it = succ_edges[ vpu_id ].find( vpv_id );
        if(edge_it == succ_edges[vpu_id].end() || edge_it->second.rce_p->weight < tmp_edge.rce_p->weight)
        {
            succ_edges[ vpu_id ][ vpv_id ] = tmp_edge;
            pred_edges[ vpv_id ][ vpu_id ] = tmp_edge;
        }
    }
}

void RevdirGraph::run_max_path_cover()
{
    vector<short> vertex_dir(succ_edges.size(), 0);

    km.n = vertex_dir.size();
    km.resize();
    km.w.assign(km.n * km.n, 0);

    INFO_PRINT(1, "run 2-approx max-cut");
    for(size_t u = 0; u < vertex_dir.size(); ++u)
    {
        long long weight_f = 0, weight_b = 0;
        for(map<size_t, GraphEdge>::iterator succ_it = succ_edges[ u ].begin();
                succ_it != succ_edges[ u ].end(); ++succ_it)
        {
            if(vertex_dir[ succ_it->first ] == 1) // succ is in R set
                weight_f += succ_it->second.rce_p->weight;
            else if(vertex_dir[ succ_it->first ] == 2) // succ is in F set
                weight_b += succ_it->second.rce_p->weight;
        }

        for(map<size_t, GraphEdge>::iterator pred_it = pred_edges[ u ].begin();
                pred_it != pred_edges[ u ].end(); ++pred_it)
        {
            if(vertex_dir[ pred_it->first ] == 1) // pred is in R set
                weight_f += pred_it->second.rce_p->weight;
            else if(vertex_dir[ pred_it->first ] == 2) // pred is in F set
                weight_b += pred_it->second.rce_p->weight;
        }

        short the_other_side = weight_f > weight_b;
        the_other_side |= !the_other_side << 1;
        vertex_dir[ u ] = the_other_side ^ 3;

        // add edge for KM
        long long uu, vv, ww;
        for(map<size_t, GraphEdge>::iterator succ_it = succ_edges[ u ].begin();
                succ_it != succ_edges[ u ].end(); ++succ_it)
        {
            if(vertex_dir[ succ_it->first ] == the_other_side)
            {
                uu = u;
                vv = succ_it->first; // the other side
                ww = succ_it->second.rce_p->weight;

                if((succ_it->second.is_type0 && the_other_side == 2) ||
                        (!succ_it->second.is_type0 && the_other_side == 1))
                {
                    swap(uu, vv);
                }
                km.w.at(uu, vv) = ww;
            }
        }

        for(map<size_t, GraphEdge>::iterator pred_it = pred_edges[ u ].begin();
                pred_it != pred_edges[ u ].end(); ++pred_it)
        {
            if(vertex_dir[ pred_it->first ] == the_other_side)
            {
                uu = pred_it->first; // the other side
                vv = u;
                ww = pred_it->second.rce_p->weight;

                if((pred_it->second.is_type0 && the_other_side == 1) ||
                        (!pred_it->second.is_type0 && the_other_side == 2))
                {
                    swap(uu, vv);
                }
                km.w.at(uu, vv) = ww;
            }
        }
    }

    DEBUG_PRINT_FILL_BLANK(1, "max-cut result:");
#ifdef DEBUG
    cerr << "    ";
    for(size_t u = 0; u < vertex_dir.size(); ++u)
        cerr << "(" << u << ", " << (vertex_dir[u] == 1 ? "R" : "F") << "); ";
    cerr << endl;
#endif

    INFO_PRINT(1, "Run KM algorithm");
    km.run();

    INFO_PRINT(1, "compute paths");
    vector<size_t> link_x( km.n );
    vector<bool> visited_x( km.n, false );
    for(size_t i = 0; i < km.n; ++i)
        link_x[ km.link_y[i] ] = i;

    for(size_t i = 0; i < link_x.size(); ++i)
    {
        if(!visited_x[i])
        {
            vector<size_t> fp, bp;
            // find forward
            bool is_cycle = false;
            size_t u = i;
            size_t min_edge_tail;

            visited_x[ u ] = true;
            while( km.w.at(u, link_x[u]) > 0)
            {
                if(fp.empty())
                    min_edge_tail = 0;
                else if(km.w.at( km.link_y[ fp[min_edge_tail] ] , fp[min_edge_tail] ) > km.w.at( u , link_x[u] ))
                    min_edge_tail = fp.size();

                u = link_x[ u ];
                fp.push_back( u );


                if(visited_x[ u ])
                {
                    is_cycle = true;
                    break;
                }

                visited_x[ u ] = true;
            }

            if(is_cycle)
            {
                paths.push_back( vector<size_t>() );
                vector<size_t>& one_path = paths.back();
                one_path.assign( fp.begin() + min_edge_tail, fp.end() );
                one_path.insert( one_path.end(),  fp.begin(), fp.begin() + min_edge_tail );
                continue;
            }

            // find backward
            u = i;
            while(km.w.at( km.link_y[u] , u) > 0)
            {
                u = km.link_y[u];
                bp.push_back( u );
                visited_x[ u ] = true;
            }

            paths.push_back( vector<size_t>() );
            vector<size_t>& one_path = paths.back();
            one_path.assign(bp.rbegin(), bp.rend());
            one_path.push_back(i);
            one_path.insert(one_path.end(), fp.begin(), fp.end());
        }
    }

    DEBUG_PRINT("output paths begin: n = %lld", km.n);
#ifdef DEBUG
    for(size_t i = 0; i < paths.size(); ++i)
    {
        cerr << "Path " << i << ": " << paths[i][0];
        for(size_t j = 1; j < paths[i].size(); ++j)
                cerr << " -> " << paths[i][j];
        cerr << endl;
    }
#endif
    DEBUG_PRINT("output paths end");
}

/*========================= class RevdirGraph =============================
================================== end ==================================*/

void reverse_nucleotide(loon::Nucleotide& nuc)
{
    swap(nuc.bases_cnt[0], nuc.bases_cnt[1]);
    swap(nuc.bases_cnt[2], nuc.bases_cnt[3]);
    switch(nuc.ch)
    {
        case 'a': case 'A':
            nuc.ch = 'T';
            break;
        case 't': case 'T':
            nuc.ch = 'A';
            break;
        case 'g': case 'G':
            nuc.ch = 'C';
            break;
        case 'c': case 'C':
            nuc.ch = 'G';
            break;
        default:
            ERROR_PRINT("Unknown bsse (%d)", int(nuc.ch));
            exit(1);
    }
}

void compute_reverse_complement(loon::Consensus& cons)
{
    size_t i = 0, j = cons.size() - 1;
    while(i < j)
    {
        swap( cons[i], cons[j] );
        reverse_nucleotide(cons[i]);
        reverse_nucleotide(cons[j]);
        ++i;
        --j;
    }
    if(i == j)
        reverse_nucleotide( cons[i] );
}

string consensus_type_reverse_string(loon::Consensus& cons)
{
    string ret;
    ret.reserve( cons.size() );
    for(loon::Consensus::const_reverse_iterator it = cons.rbegin(); it != cons.rend(); ++it)
    {
        switch(it->ch)
        {
            case 'a': case 'A':
                ret.push_back('T');
                break;
            case 't': case 'T':
                ret.push_back('A');
                break;
            case 'g': case 'G':
                ret.push_back('C');
                break;
            case 'c': case 'C':
                ret.push_back('G');
                break;
            default:
                ERROR_PRINT("Unknown bsse (%d)", int(it->ch));
                exit(1);
        }
    }
    return ret;
}

void do_consensus_save_path_head(bool is_type0, 
                                Consensus& ctg,
                                vector<Consensus>& revdir_res)
{
    DEBUG_PRINT("add path head");
    if(!is_type0)
    {
        compute_reverse_complement( ctg );
        ctg.reverse_path();
    }
    revdir_res.push_back( Consensus() );
    revdir_res.back().swap( ctg );
}

void do_consensus_save_edge(bool keep_edge_dir,
                            bool is_type0, 
                            loon::Consensus& edge, 
                            Consensus& succ,
                            vector<Consensus>& revdir_res)
{
    // consensus edge
    if(keep_edge_dir)
    {
        try
        {
            DEBUG_PRINT("consensus [1]");
            if(run_options & 0x10)
                revdir_res.back().consensus_right( loon::blasr_align_choice, edge.toString() );
            else
                revdir_res.back().consensus_right( edge.toString() );
        }catch(int error_code)
        {
            DEBUG_PRINT("consensus error [1]: error_code = %d", error_code);
            revdir_res.push_back( edge );
        }
    }
    else
    {
        try
        {
            DEBUG_PRINT("consensus [2]");
            if(run_options & 0x10)
                revdir_res.back().consensus_right( loon::blasr_align_choice, consensus_type_reverse_string( edge ) );
            else
                revdir_res.back().consensus_right( consensus_type_reverse_string( edge ) );
        }catch(int error_code)
        {
            DEBUG_PRINT("consensus error [2]: error_code = %d", error_code);
            revdir_res.push_back( edge );
            compute_reverse_complement( revdir_res.back() );
        }
    }

    // consensus succ
    if(is_type0)
    {
        compute_reverse_complement( succ );
        succ.reverse_path();
    }
    try
    {
        DEBUG_PRINT("consensus [3]");
        if(run_options & 0x10)
            revdir_res.back().consensus_right( loon::blasr_align_choice, succ.toString(),
                revdir_res.back().size() - min(edge.size(), succ.size()));
        else
            revdir_res.back().consensus_right( succ.toString(),
                revdir_res.back().size() - min(edge.size(), succ.size()));
        revdir_res.back().append_forward( succ );
    }catch(int error_code)
    {
        DEBUG_PRINT("consensus error [3]: error_code = %d", error_code);
        revdir_res.push_back( loon::Consensus() );
        revdir_res.back().swap( succ );
    }
}

void bridge_revdir(vector<Consensus>& contigs,
        const loon::Multiseq& long_reads,
        vector<Consensus>& revdir_res,
        vector<bool>& used_long_reads)
{
    const size_t lr_n = long_reads.size();
    used_long_reads.assign(lr_n, false);
    RPrimalGraph rp_graph;

    loon::Aligner aligner;
    aligner.set_simple_align();

    INFO_PRINT("Convert contigs to strings");
    loon::Multiseq tmp_ref;
    tmp_ref.resize(contigs.size());
    for(size_t i = 0; i < contigs.size(); ++i)
        tmp_ref[i].seq = contigs[i].toString();

    INFO_PRINT("All contigs vs. all long reads matching");
    nucmer_Deltas deltas_forward, deltas_backward;
    loon::all_v_all_matching(tmp_ref, long_reads, deltas_forward, deltas_backward, true);

    vector<vector<size_t> > fleft_vertices( lr_n ), fright_vertices( lr_n );
    vector<vector<size_t> > rleft_vertices( lr_n ), rright_vertices( lr_n );

    // compute vertices for edges
    INFO_PRINT("Compute vertices for edges ");
    progress.start( deltas_forward.size() + deltas_backward.size() );
    for(size_t ii = 0; ii < deltas_forward.size(); ++ii)
    {
        size_t lr_i = str2lu( deltas_forward[ ii ].qry_tag );
        used_long_reads[lr_i] = add_edge(fleft_vertices[ lr_i ], 
                                            fright_vertices[ lr_i ], 
                                            deltas_forward[ ii ]) || 
                                    used_long_reads[lr_i];

        progress.progress(1);
    }
    for(size_t ii = 0; ii < deltas_backward.size(); ++ii)
    {
        size_t lr_i = str2lu( deltas_backward[ ii ].qry_tag );
        used_long_reads[lr_i] = add_edge(rleft_vertices[ lr_i ], 
                                            rright_vertices[ lr_i ], 
                                            deltas_backward[ ii ]) || 
                                    used_long_reads[lr_i];

        progress.progress(1);
    }
    progress.stop();

    // build the revdir primal graph
    INFO_PRINT("Build the revdir primal graph");

    progress.start( lr_n );
    for(size_t lr_i = 0; lr_i < lr_n; ++lr_i)
    {
        if(!fleft_vertices[ lr_i ].empty() && !rleft_vertices[ lr_i ].empty())
        {
            rp_graph.add_type0_edge(long_reads[ lr_i ], 
                                    long_reads[ lr_i ].seq, 
                                    long_reads[ lr_i ].rev, 
                                    fleft_vertices[ lr_i ], 
                                    rleft_vertices[ lr_i ]);
        }
        
        if(!fright_vertices[ lr_i ].empty() && !rright_vertices[ lr_i ].empty())
        {
            rp_graph.add_type1_edge(long_reads[ lr_i ], 
                                    long_reads[ lr_i ].seq,
                                    long_reads[ lr_i ].rev, 
                                    fright_vertices[ lr_i ], 
                                    rright_vertices[ lr_i ]);
        }

        progress.progress( 1 );
    }

    DEBUG_PRINT("Print primal graph");
#ifdef DEBUG
    for(map<VertexPair, set<size_t> >::iterator me_it = rp_graph.multiedges.begin();
            me_it != rp_graph.multiedges.end(); ++me_it)
    {
        DEBUG_PRINT("vertex (u, v) = (%lu, %lu), edge_type = %d, edges (id, weight): ", 
                me_it->first.u, me_it->first.v, !(me_it->first.is_type0));
        for(set<size_t>::iterator e_it = me_it->second.begin(); e_it != me_it->second.end(); ++e_it)
        {
            cerr << "(" << (*e_it) << ", " << rp_graph.rcom_edges[ *e_it ].weight << "); ";
        }
        cerr << endl;
    }
#endif

    INFO_PRINT("Build the actual graph based on the revdir primal graph");
    RevdirGraph revdir_graph(contigs.size(), rp_graph, contigs);

    DEBUG_PRINT("Print actual graph");
    DEBUG_PRINT("Vertex translation: fx to x");
#ifdef DEBUG
    for(size_t fx = 0; fx < revdir_graph.vertex_fxtox.size(); ++fx)
        DEBUG_PRINT("fx = %lu\tx = %lu", fx, revdir_graph.get_vertex_original_id(fx));
    DEBUG_PRINT("The graph");
    for(size_t x = 0; x < revdir_graph.succ_edges.size(); ++x)
    {
        DEBUG_PRINT("vertex [%lu], succ = (vertex id, weight, edge type)", x);
        for(map<size_t, RevdirGraph::GraphEdge>::iterator succ_it = revdir_graph.succ_edges[x].begin();
                succ_it != revdir_graph.succ_edges[x].end(); ++succ_it)
        {
            cerr << "(" << (succ_it->first) << ", " << (succ_it->second.rce_p->weight) << ", " << (!(succ_it->second.is_type0)) << "); ";
        }
        cerr << endl;
    }
#endif

    INFO_PRINT("Run max cut & max path cover (KM) algorithm");
    revdir_graph.run_max_path_cover();

    INFO_PRINT("Do consensus");
    size_t fx_uid, fx_vid;
    RComEdge* rce_p;
    bool is_type0;
    bool keep_edge_dir;
    for(size_t i = 0; i < revdir_graph.paths.size(); ++i)
    {
        fx_uid = revdir_graph.paths[i][0];
        if(revdir_graph.paths[i].size() == 1)
        {
            DEBUG_PRINT("single-vertex path");
            revdir_res.push_back( Consensus() );
            revdir_res.back().swap( contigs[ revdir_graph.get_vertex_original_id( fx_uid ) ] );
        #ifdef DEBUG
            cerr << endl;
        #endif
        }
        else
        {
            fx_vid = revdir_graph.paths[i][1];
            try
            {
                keep_edge_dir = revdir_graph.find_edge(rce_p, is_type0, fx_uid, fx_vid);
            }catch(int error_code)
            {
                ERROR_PRINT("cannot find edge");
                exit(1);
            }
            
            do_consensus_save_path_head(is_type0, contigs[ revdir_graph.get_vertex_original_id( fx_uid ) ],
                    revdir_res);
            do_consensus_save_edge(keep_edge_dir, is_type0,
                    rce_p->rcom_edge, contigs[ revdir_graph.get_vertex_original_id( fx_vid ) ],
                    revdir_res);
        #ifdef DEBUG
            cerr << endl;
        #endif

            for(size_t j = 2; j < revdir_graph.paths[i].size(); ++j)
            {
                fx_uid = fx_vid;
                fx_vid = revdir_graph.paths[i][j];
                try
                {
                    keep_edge_dir = revdir_graph.find_edge( rce_p, is_type0, fx_uid, fx_vid );
                }catch(int error_code)
                {
                    ERROR_PRINT("cannot find edge");
                    exit(1);
                }
                do_consensus_save_edge(keep_edge_dir, is_type0,
                        rce_p->rcom_edge, contigs[ revdir_graph.get_vertex_original_id( fx_vid ) ],
                        revdir_res);
            #ifdef DEBUG
                cerr << endl;
            #endif
            }
        }
    }

    INFO_PRINT("Collect unused vertices");
    for(size_t i = 0; i < revdir_graph.n; ++i)
    {
        if(revdir_graph.vertex_xtofx[ i ] == loon::INF)
        {
            revdir_res.push_back( Consensus() );
            revdir_res.back().swap( contigs[ i ] );
        }
    }
}

template<class T>
void keep_chosen_vector(T& a, const vector<size_t>& kept_idx)
{
    size_t a_id = 0;
    for(vector<size_t>::const_iterator cit = kept_idx.begin(); cit != kept_idx.end(); ++cit)
    {
        if(a_id < *cit)
            a[ a_id ] = a[ *cit];
        ++a_id;
    }
    a.resize(a_id);
}

template<class T>
void keep_chosen_vector(T& a, const vector<bool>& kept_idx, const bool keep_opt = true)
{
    size_t a_id = 0;
    for(size_t i = 0; i < kept_idx.size(); ++i)
        if(kept_idx[ i ] == keep_opt)
        {
            if(a_id < i)    a[a_id] = a[i];
            ++a_id;
        }
    a.resize( a_id );
}

// void do_bridging_debug(vector<loon::Consensus>& contigs, loon::Multiseq& qry, vector<size_t>& lr_idx)
// {
//     if(contigs.size() <= 1)
//         return;
// 
//     vector<loon::Consensus> ctg_res;
//     vector<bool> lr_used;
// 
//     keep_chosen_vector(qry, lr_idx);
//     size_t round_cnt = 0;
//     {
//         INFO_PRINT("Bridge co-direction contigs, round #%lu", ++round_cnt);
//         bridge_codir(contigs, qry, ctg_res, lr_used);
// 
//         contigs.swap( ctg_res );
//         if(contigs.size() <= 1)
//             return;
//         size_t last_qry_size = qry.size();
//         keep_chosen_vector(qry, lr_used);
//         if(last_qry_size == qry.size())
//             return;
//         ctg_res.clear();
//         INFO_PRINT("bridge co-direction contigs done! exit");
// 
//         // INFO_PRINT("Bridge rev-direction contigs, round #%lu", round_cnt);
//         // bridge_revdir(contigs, qry, ctg_res, lr_used);
// 
//         // contigs.swap( ctg_res );
//         // if(contigs.size() <= 1)
//         //     return;
//         // last_qry_size = qry.size();
//         // keep_chosen_vector(qry, lr_used);
//         // if(last_qry_size == qry.size())
//         //     return;
//         // ctg_res.clear();
//     }
// }

void do_bridging(vector<Consensus>& contigs, loon::Multiseq& qry, vector<size_t>& lr_idx)
{
    if(contigs.size() <= 1)
        return;

    vector<Consensus> ctg_res;
    vector<bool> lr_used;

    keep_chosen_vector(qry, lr_idx);
    size_t round_cnt = 0;
    while(true)
    {
        INFO_PRINT("Bridge co-direction contigs, round #%lu", ++round_cnt);
        bridge_codir(contigs, qry, ctg_res, lr_used);

        contigs.swap( ctg_res );
        if(contigs.size() <= 1)
            return;
        size_t last_qry_size = qry.size();
        keep_chosen_vector(qry, lr_used);
        if(last_qry_size == qry.size())
            return;
        ctg_res.clear();

        INFO_PRINT("Bridge rev-direction contigs, round #%lu", round_cnt);
        bridge_revdir(contigs, qry, ctg_res, lr_used);

        contigs.swap( ctg_res );
        if(contigs.size() <= 1)
            return;
        last_qry_size = qry.size();
        keep_chosen_vector(qry, lr_used);
        if(last_qry_size == qry.size())
            return;
        ctg_res.clear();
    }
}

/*================== class NCAlignment ==================
========================= start =======================*/
class NCAlignment
{
public:
    size_t ref_start, ref_end;
    string* strp;
public:
    NCAlignment(size_t rs = 0, size_t rt = 0, string* sp = NULL):
        ref_start(rs), ref_end(rt), strp(sp)
    {}
};
/*================== class NCAlignment ==================
========================= end ==========================*/

/*================== class Region =====================
======================= start =======================*/
class Region
{
public:
    size_t start, end;
    list<NCAlignment>* long_reads;
public:
    Region();
    void sort();
    static bool cmp_pos(const NCAlignment &nc1, const NCAlignment &nc2);
};

Region::Region(): start(0), end(0), long_reads(NULL)
{
}

void Region::sort()
{
    if(long_reads)
        long_reads->sort( cmp_pos );
}

/*static*/ bool Region::cmp_pos(const NCAlignment &nc1, const NCAlignment &nc2)
{
    return (nc1.ref_start < nc2.ref_start ||
            (nc1.ref_start == nc2.ref_start && 
             nc1.ref_end > nc2.ref_end));
}
/*================== class Region =====================
======================== end ========================*/

/*================== class Regions ====================
======================== start ======================*/
class Regions: public vector<Region>
{
public:
    long long find_interval(size_t st, size_t ed, long long start_it = 0);
    void add_coverage(vector<tigrinc::Nucmer_Cluster>& c, string* strp);
    void sort();
};

long long Regions::find_interval(size_t st, size_t ed, long long start_it/* = 0*/)
{
    for(; start_it < this->size(); ++start_it)
    {
        if(this->at( start_it ).start + MIN_OVERLAP > ed ||
                this->at( start_it ).end < st + MIN_OVERLAP)
            continue;
        return start_it;
    }
    return -1;
}

void Regions::add_coverage(vector<tigrinc::Nucmer_Cluster>& c, string* strp)
{
    if(c.empty())
        return;

    Region res;
    res.start = c.front().ref_start;
    res.end = c.back().ref_end;
    res.long_reads = new list<NCAlignment>();
    res.long_reads->push_back( NCAlignment(res.start, res.end, strp) );

    long long pos = find_interval(res.start, res.end);
    while(pos != -1)
    {
        res.start = min( res.start, this->at(pos).start);
        res.end = max( res.end, this->at(pos).end );
        res.long_reads->splice( res.long_reads->end(), *(this->at(pos).long_reads) );
        delete this->at(pos).long_reads;
        this->at(pos) = this->back();
        this->pop_back();

        pos = find_interval(res.start, res.end, pos);
    }

    this->push_back(res);

}

void Regions::sort()
{
    for(vector<Region>::iterator rgn_it = this->begin(); rgn_it != this->end(); ++rgn_it)
        rgn_it->sort();
}
/*================== class Regions ====================
======================== end ========================*/

double compute_delta_score(const vector<tigrinc::Nucmer_Cluster>& clu, const long long len)
{
    const size_t max_diff = 30;

    if(clu.empty())
        return 0;

    long long cov = 0;
    long long sep = 0;
    long long last_ref_end = clu[0].ref_start - 1;
    long long last_qry_end = clu[0].qry_start - 1;
    long long ref_gap, qry_gap;
    for(size_t i = 0; i < clu.size(); ++i)
    {
        ref_gap = clu[i].ref_start - last_ref_end - 1;
        qry_gap = clu[i].qry_start - last_qry_end - 1;
        if(abs(ref_gap - qry_gap) > max_diff)
            return 0;
        
        cov += clu[i].qry_end - clu[i].qry_start + 1;
        sep += ref_gap + qry_gap;

        last_ref_end = clu[i].ref_end;
        last_qry_end = clu[i].qry_end;
    }

    return alpha * cov / len + beta * len / (sep + len + len);
}

void open_file(const char* file_suffix, ofstream& fout, bool binary = false)
{
    string fname = prefix + file_suffix;
    if(binary)
        fout.open( fname.c_str(), ofstream::binary);
    else
        fout.open( fname.c_str() );
    loon::check_file_open( fout, fname );
}

void read_contigs(vector<Consensus>& contigs)
{
    ifstream fin(ctg_file.c_str(), ifstream::binary);
    loon::check_file_open( fin, ctg_file );

    size_t cnt = 0;
    int c;
    while((c = fin.peek()) != EOF)
    {
        contigs.push_back( Consensus() );
        contigs.back().read_binary( fin );
        contigs.back().append_segment( cnt );
        ++cnt;
    }

    fin.close();

    string segfile = prefix + ".dis";
    fin.open(segfile.c_str());
    loon::check_file_open(fin, segfile);

    SegInfo tmp_seginfo;
    while(fin >> tmp_seginfo.ref_id >> tmp_seginfo.start_loc >> tmp_seginfo.end_loc)
        seginfo_vec.push_back(tmp_seginfo);

    if(seginfo_vec.size() != cnt)
    {
        ERROR_PRINT("[%s] and [%s] don\'t match", segfile.c_str(), ctg_file.c_str());
        exit(1);
    }

    fin.close();
}

void read_idx(vector<size_t>& lr_idx)
{
    ifstream fin(idx_file.c_str());
    loon::check_file_open( fin, idx_file );

    size_t n;
    fin >> n;
    lr_idx.resize(n);

    for(size_t i = 0; i < n; ++i)
        fin >> lr_idx[ i ];

    fin.close();
}

void save_final_result(vector<Consensus>& contigs)
{
    ofstream fout_finalbin, fout_finalfasta, fout_path;
    open_file( ".final.bin", fout_finalbin, true );
    open_file( ".final.fasta", fout_finalfasta);
    open_file( ".path", fout_path);

    loon::Fasta tmp_fasta;
    for(size_t i = 0; i < contigs.size(); ++i)
    {
        contigs[i].save_binary( fout_finalbin );

        ostringstream oss;
        oss << "cons_" << i;
        tmp_fasta.tag = oss.str();
        tmp_fasta.seq = contigs[i].toString();
        tmp_fasta.save_fasta( fout_finalfasta );

        contigs[ i ].print_path( fout_path, i );
    }
    
    fout_finalbin.close();
    fout_finalfasta.close();
    fout_path.close();

    INFO_PRINT("%lu contigs", contigs.size());

#ifdef DEBUG_CHGCNT
    INFO_PRINT("%lu changes", chgcnt - contigs.size());
#endif
}

bool deltas_tag_cmp(const tigrinc::Nucmer_Delta& d1, const tigrinc::Nucmer_Delta& d2)
{
    int cmp1 = d1.ref_tag.compare( d2.ref_tag );
    return (cmp1 < 0 || (cmp1 == 0 && d1.qry_tag < d2.qry_tag));
}


void try_add_regions(tigrinc::Nucmer_Delta* df, tigrinc::Nucmer_Delta *db, 
        loon::Multiseq& qry, Regions& regions,
        vector<bool>& ref_used, vector<bool>& qry_used)
{
    double f_score = ( df == NULL ? 0 : compute_delta_score(df->clusters, df->qry_len));
    double b_score = ( db == NULL ? 0 : compute_delta_score(db->clusters, db->qry_len));
    if(f_score >= b_score)
    {
        if(f_score >= min_score)
        {
            size_t qry_id = str2lu( df->qry_tag );
            regions.add_coverage( df->clusters, &qry[ qry_id ].seq);
            qry_used[ qry_id ] = true;
            ref_used[ str2lu(df->ref_tag) ] = true;
        }
    }
    else
    {
        if(b_score >= min_score)
        {
            size_t qry_id = str2lu( db->qry_tag );
            regions.add_coverage(db->clusters, &qry[ qry_id ].rev);
            qry_used[ qry_id ] = true;
            ref_used[ str2lu(db->ref_tag) ] = true;
        }
    }
}

void consensus_regions(Regions& regions, vector<Consensus>& contigs,
        ostream& fout_cons_fasta, ostream& fout_cons_bin, const string& ref_seq, size_t ref_id)
{
    loon::Fasta tmp_fasta;
    if(regions.empty())
        return;

    INFO_PRINT("Consensus regions for a reference: %lu regions", regions.size());
    loon::Progress progress;
    progress.start(regions.size());
    
    size_t seg_id, start_loc, end_loc;
    for(Regions::iterator rgn_it = regions.begin(); rgn_it != regions.end(); ++rgn_it)
    {
        rgn_it->sort();
        // three strategies
        if(consensus_strategy == 0)
        {
            seg_id = seginfo_vec.size();
            start_loc = rgn_it->start;
            end_loc = rgn_it->end;
            seginfo_vec.push_back( SegInfo(ref_id, start_loc, end_loc) );

            contigs.push_back( Consensus() );
            contigs.back().set_str( ref_seq.substr(start_loc - 1, end_loc - start_loc + 1) );
            contigs.back().append_segment(seg_id);
        }
        else
        {
            list<NCAlignment>::iterator last_nca;

            contigs.push_back( Consensus() );
            Consensus& cons = contigs.back();

            for(list<NCAlignment>::iterator nc_it = rgn_it -> long_reads -> begin();
                    nc_it !=  rgn_it -> long_reads -> end(); ++nc_it)
            {
                if(cons.empty())
                {
                    seg_id = seginfo_vec.size();
                    start_loc = nc_it->ref_start;
                    end_loc = nc_it->ref_end;
                    seginfo_vec.push_back( SegInfo(ref_id, start_loc, end_loc) );

                    cons.set_str( *nc_it->strp );
                    if(consensus_strategy == 1)
                        last_nca = nc_it;
                    cons.append_segment( seg_id );
                }
                else
                {
                    try
                    {
                        if(consensus_strategy == 1)
                        {
                            if(last_nca -> ref_end >= nc_it -> ref_end)
                                continue;
                            last_nca = nc_it;
                        }
                        if(run_options & 0x10)
                            cons.consensus_right_auto_delimiter(loon::blasr_align_choice, *nc_it->strp, false );
                        else
                            cons.consensus_right_auto_delimiter( *nc_it->strp, false );

                        seginfo_vec.back().set_end_loc( nc_it->ref_end );

                    }catch(int e)
                    {
                        fout_log << "Consensus error[2]: error_code = " << e << endl;
                    }
                }
            }
        }

        if(run_options & 1)
        {
            Consensus& cons = contigs.back();
            cons.save_binary( fout_cons_bin );

            tmp_fasta.seq = cons.toString();
            ostringstream oss;
            oss << "cons_" << contigs.size();
            tmp_fasta.tag = oss.str();
            tmp_fasta.save_fasta( fout_cons_fasta );
        }
    #ifdef DEBUG_CHGCNT
        ++chgcnt;
    #endif

        progress.progress(1);
    }
    regions.clear();
    progress.stop();
}

void do_alignment_and_consensus(loon::Multiseq& ref, loon::Multiseq& qry,
        vector<Consensus>& contigs, vector<size_t>& lr_idx)
{
    size_t should_have_regions_cnt = 0;
    ofstream fout_cons_fasta, fout_cons_bin, fout_idx;
    if(run_options & 1)
    {
        open_file(".cons.fasta", fout_cons_fasta);
        open_file(".cons.bin", fout_cons_bin, true);
        open_file(".idx.txt", fout_idx);
    }

    INFO_PRINT("Align all the long reads to the reference genome");
    nucmer_Deltas deltas_forward, deltas_backward;
    if(filter_prefix.empty())
    {
        loon::all_v_all_matching(ref, qry, deltas_forward, deltas_backward);
    }
    else
    {
        loon::read_deltas(deltas_forward, filter_prefix + ".forward.filter");
        loon::read_deltas(deltas_backward, filter_prefix + ".backward.filter");
        deltas_forward.trim();
        deltas_backward.trim();
    }
    sort(deltas_forward.begin(), deltas_forward.end(), deltas_tag_cmp);
    sort(deltas_backward.begin(), deltas_backward.end(), deltas_tag_cmp);

    Regions regions;
    size_t deltas_fn = deltas_forward.size();
    size_t deltas_bn = deltas_backward.size();
    size_t i=0 , j=0;
    string last_tag = "";
    vector<bool> qry_used( qry.size(), false );
    vector<bool> ref_used( ref.size(), false );

    INFO_PRINT("Compute intervals");
    progress.start(deltas_fn + deltas_bn);
    while(i < deltas_fn && j < deltas_bn)
    {
        if(deltas_forward[i].ref_tag == last_tag)
        {
            if(deltas_backward[j].ref_tag == last_tag)
            {// ftag == last_tag == btag
                int cmp_tag2 = deltas_forward[i].qry_tag.compare( deltas_backward[j].qry_tag );
                if(cmp_tag2 < 0)
                {
                    try_add_regions(&deltas_forward[i], NULL, qry, regions, ref_used, qry_used);
                    ++i;

                    progress.progress(1);
                }
                else if(cmp_tag2 > 0)
                {
                    try_add_regions(NULL, &deltas_backward[j], qry, regions, ref_used, qry_used);
                    ++j;
                    progress.progress(1);
                }
                else
                {
                    try_add_regions(&deltas_forward[i], &deltas_backward[j], qry, regions, ref_used, qry_used);
                    ++i;    ++j;
                    progress.progress(2);
                }
            }
            else
            {// ftag == last_tag && last_tag < btag
                try_add_regions(&deltas_forward[i], NULL, qry, regions, ref_used, qry_used);
                ++i;
                progress.progress(1);
            }
        }
        else if(deltas_backward[j].ref_tag == last_tag)
        {// ftag > last_tag && btag == last_tag
            try_add_regions(NULL, &deltas_backward[j], qry, regions, ref_used, qry_used);
            ++j;
            progress.progress(1);
        }
        else
        {// ftag > last_tag && btag > last_tag
            should_have_regions_cnt += regions.size();
            size_t ref_id = str2lu(last_tag);
            consensus_regions(regions, contigs, fout_cons_fasta, fout_cons_bin, ref[ ref_id ].seq, ref_id);
            last_tag = min(deltas_forward[i].ref_tag, deltas_backward[j].ref_tag);
        }
    }

    while(i < deltas_fn)
    {
        if(deltas_forward[i].ref_tag != last_tag)
        {
            should_have_regions_cnt += regions.size();
            size_t ref_id = str2lu( last_tag );
            consensus_regions(regions, contigs, fout_cons_fasta, fout_cons_bin, ref[ ref_id ].seq, ref_id);
            last_tag = deltas_forward[i].ref_tag;
        }

        try_add_regions(&deltas_forward[i], NULL, qry, regions, ref_used, qry_used);
        ++i;
        progress.progress(1);
    }

    while(j < deltas_bn)
    {
        if(deltas_backward[j].ref_tag != last_tag)
        {
            should_have_regions_cnt += regions.size();
            size_t ref_id = str2lu(last_tag);
            consensus_regions(regions, contigs, fout_cons_fasta, fout_cons_bin, ref[ ref_id ].seq, ref_id);
            last_tag = deltas_backward[j].ref_tag;
        }
        try_add_regions(NULL, &deltas_backward[j], qry, regions, ref_used, qry_used);
        ++j;
        progress.progress(1);
    }
    progress.stop();
    should_have_regions_cnt += regions.size();
    size_t ref_id = str2lu(last_tag);
    consensus_regions(regions, contigs, fout_cons_fasta, fout_cons_bin, ref[ ref_id ].seq, ref_id);
    INFO_PRINT("There should be [%lu] regions", should_have_regions_cnt);

    if(run_options & 8)
    {// don't bridge the unaligned references
        string fname = prefix + ".una.fa";
        ofstream fout_unaligned(fname.c_str());
        loon::check_file_open(fout_unaligned, fname);

        for(size_t i = 0; i < ref.size(); ++i)
            if(!ref_used[ i ])
                ref[i].save_fasta( fout_unaligned );

        fout_unaligned.close();
    }
    else
    {// bridge the unaligned references too
        for(size_t i = 0; i < ref.size(); ++i)
        {
            if(!ref_used[ i ])
            {
                seginfo_vec.push_back( SegInfo(i, 1, ref[i].seq.length()) );

                contigs.push_back( Consensus() );
                contigs.back().set_str( ref[i].seq );
                contigs.back().append_segment( seginfo_vec.size() );

                if(run_options & 1)
                {
                    contigs.back().save_binary( fout_cons_bin );
                    ref[i].save_fasta( fout_cons_fasta );
                }
            }
        }
    }

    // save unused qry
    for(size_t j = 0; j < qry_used.size(); ++j)
    {
        if(!qry_used[j])
            lr_idx.push_back( j );
    }

    if(run_options & 1)
    {
        // save idx
        fout_idx << lr_idx.size() << endl;
        for(vector<size_t>::iterator it = lr_idx.begin(); it != lr_idx.end(); ++it)
            fout_idx << (*it) << ' ';
        fout_idx << endl;

        // close files
        fout_cons_fasta.close();
        fout_cons_bin.close();
        fout_idx.close();
    }

    // save disassembling info
    ofstream fout_seginfo;
    open_file(".dis", fout_seginfo);
    for(size_t i = 0; i < seginfo_vec.size(); ++i)
    {
        fout_seginfo << seginfo_vec[i].ref_id << ' ' 
                    << seginfo_vec[i].start_loc << ' ' 
                    << seginfo_vec[i].end_loc << endl;
    }
    fout_seginfo.close();
}

void help_msg(char* exec)
{
    cerr << "Usage: " << exec << " [options] <ref_file> <qry_file>" << endl;
    cerr << "Options:" << endl;
    cerr << "\t-p prefix        set prefix for output files (default smsc)" << endl;
    cerr << "\t-E               choose Blasr as alignment engine (default nucmer)" << endl;
    cerr << "\t-m               do alignment and consensus step" << endl;
    cerr << "\t-b               do bridge step" << endl;
    cerr << "\t-c congtig_file  contig file for bridge step" << endl;
    cerr << "\t-i index_file    long reads index file for bridge step" << endl;
    cerr << "\t-f filter_prefix prefix for delta-filter file" << endl;
    cerr << "\t-C strategy_id   choose consensus (default 1)" << endl;
    cerr << "\t                     * 0: don\'t consensus, use the blocks in the ref directly" << endl;
    cerr << "\t                     * 1: do simple consensus" << endl;
    cerr << "\t                     * 2: do complicated consensus (result is bad ToT)" << endl;
    cerr << "\t-B               Don't bridge the unaligned references" << endl;
    cerr << "\t-j threads       # of threads" << endl;
    cerr << "\t-h               show help info" << endl;
    cerr << "Notice:" << endl;
    cerr << "\tThe -c and -i works only when only do bridging" << endl;
}

void parse_parameters(int argc, char* argv[])
{
    int opt;
    istringstream iss;
    while((opt = getopt(argc, argv, "p:Embc:C:Bi:f:j:h")) != -1)
    {
        switch(opt)
        {
            case 'p':
                prefix = optarg;
                break;
            case 'E':
                run_options |= 0x10;
                break;
            case 'm':
                run_options |= 1;
                break;
            case 'b':
                run_options |= 2;
                break;
            case 'c':
                ctg_file = optarg;
                break;
            case 'C':
                consensus_strategy = optarg[0] - '0';
                if(consensus_strategy < 0 || consensus_strategy > 2)
                {
                    ERROR_PRINT("strategy_id should be 0, 1, or 2");
                    exit(1);
                }
                break;
            case 'B':
                run_options |= 8;
                break;
            case 'i':
                idx_file = optarg;
                break;
            case 'f':
                filter_prefix = optarg;
                break;
            case 'j':
                n_cores = str2lu( optarg );
                break;
            case 'h':
                help_msg( argv[0] );
                exit(0);
            default:
                ERROR_PRINT("Error: Unknown options: %c", char(opt));
                help_msg( argv[0] );
                exit(1);
        }
    }

    if(optind >= argc)
    {
        ERROR_PRINT("Error: <ref_file> and <qry_file> are not specified");
        help_msg( argv[0] );
        exit(1);
    }

    ref_file = argv[ optind++ ];

    if(optind >= argc)
    {
        ERROR_PRINT("Error: <qry_file> is not specified");
        help_msg( argv[0] );
        exit(1);
    }

    qry_file = argv[ optind++ ];
    if(optind < argc)
    {
        ERROR_PRINT("Error: Too many parameters");
        help_msg( argv[0] );
        exit(1);
    }


    if(run_options & 3)
        run_options &= 0xfffb;

    if((run_options == 2) && (ctg_file.empty() || idx_file.empty()))
    {
        ERROR_PRINT("-c and -i options should be used when only do bridging");
        help_msg( argv[0] );
        exit(1);
    }

}

string build_tmp_path()
{
    ostringstream oss;
    char *path = (char *)malloc(256), *ret;
    size_t size = 256;
    if(path == NULL)
    {
        ERROR_PRINT("Allocate memory for current path failed!");
    }

    while( (ret = getcwd(path, size)) == NULL)
    {
        if(errno == ERANGE)
        {
            size <<= 1;
            if((path = (char *)realloc(path, size)) == NULL)
            {
                ERROR_PRINT("Reallocate memory for current path failed!");
            }
        }
        else
        {
            ERROR_PRINT("get current path error!!!");
        }
    }
    oss << path << "/tmp_mummer_api_" << getpid();
    free(path);
    return oss.str();
}

int main(int argc, char* argv[])
{
    parse_parameters(argc, argv);

    DEBUG_PRINT("ref_file = [%s]", ref_file.c_str());
    DEBUG_PRINT("qry_file = [%s]", qry_file.c_str());
    DEBUG_PRINT("prefix = [%s]", prefix.c_str());
    DEBUG_PRINT("ctg_file = [%s]", ctg_file.c_str());
    DEBUG_PRINT("idx_file = [%s]", idx_file.c_str());
    DEBUG_PRINT("run_options = %d", run_options);
    DEBUG_PRINT("maximum query length is %lu", (Showuint) ~((Uint)0));

    loon::mumapi_set_script_path();
    loon::mumapi_setpoolsize(n_cores);
    loon::mumapi_setdir(build_tmp_path());

    //loon::mumapi_setdir("/Users/Loon/loon/debug/rgeca/v5.0/gambiae/tmp_mummer_api_19990");

    loon::mumapi_mkdir();

    open_file(".log", fout_log);
    INFO_PRINT("read ref_file and qry_file");
    loon::Multiseq ref, qry;
    ref.read_fasta( ref_file.c_str() );
    qry.read_fasta( qry_file.c_str() );

    INFO_PRINT("# ref = %lu, # qry = %lu", ref.size(), qry.size());

    INFO_PRINT("compute reverse complement of qry");
    qry.compute_reverse_complement();

    vector<Consensus> contigs;
    vector<size_t> lr_idx;

    if(run_options & 0x10)
    {
        INFO_PRINT("Use [Blasr] as alignment engine");
        loon::mumapi_set_engine( BLASR_ENGINE );
    }
    else
    {
        INFO_PRINT("Use [Nucmer] as alignment engine");
    }

    if(run_options & 5)
    {
        if(run_options & 1)
        {
            INFO_PRINT("do alignment and consensus; and save the result");
        }
        else
        {
            INFO_PRINT("do alignment and consensus");
        }
        do_alignment_and_consensus(ref, qry, contigs, lr_idx);
    }
    
    if(run_options & 6)
    { // do bridging step
        if(run_options == 2)
        { // read contigs file and index file
            INFO_PRINT("read contigs file and index file");

            read_contigs(contigs);
            read_idx(lr_idx);
        }

        INFO_PRINT("do bridging");
        do_bridging(contigs, qry, lr_idx);

        INFO_PRINT("save bridging result");
        save_final_result(contigs);
    }

    loon::mumapi_rmdir();

    // INFO_PRINT("read contigs file and index file");
    // read_contigs(contigs);
    // read_idx(lr_idx);
    // INFO_PRINT("do bridging");
    // do_bridging_debug(contigs, qry, lr_idx);


    INFO_PRINT("done!");
    fout_log.close();
    return 0;
}
