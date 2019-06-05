#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <iterator>
#include <list>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <errno.h>

#include "path_cover.h"

#include "common-bio/Progress.h"
#include "common-bio/Consensus.h"
#include "common-bio/Aligner.h"
#include "common-bio/Multiseq.h"
#include "common-bio/util.h"
#include "common-bio/baseheader.h"


#include "nucmer.h"

#include "mummer_api.h"
#include "msa_consensus.h"

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

inline size_t head_idx(size_t s)
{
    return (s << 1);
}

inline size_t tail_idx(size_t s)
{
    return (s << 1 | 1);   
}

inline size_t forward_eid(size_t e)
{
    return (e << 1);
}

inline size_t reverse_eid(size_t e)
{
    return (e << 1 | 1);
}

size_t str2lu(const string& str)
{
    size_t ret = 0;
    for(size_t i = 0; i < str.length(); ++i)
        ret = ret * 10 + str[i] - '0';
    return ret;
}

bool deltas_tag_cmp(const tigrinc::Nucmer_Delta& d1, const tigrinc::Nucmer_Delta& d2)
{
    int cmp1 = d1.ref_tag.compare( d2.ref_tag );
    return (cmp1 < 0 || (cmp1 == 0 && d1.qry_tag < d2.qry_tag));
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
// class SegInfo
// {
// public:
//     size_t ref_id, start_loc, end_loc;
// public:
//     SegInfo(size_t r_id = 0, size_t s_loc = 1, size_t e_loc = 0):
//         ref_id(r_id), start_loc(s_loc), end_loc(e_loc)
//     {
//         if(s_loc == 0)
//         {
//             ERROR_PRINT("Implementation Error: start_loc starts from 0, not from 1");
//             exit(1);
//         }
//     }
// 
//     void set_ref_id(size_t r_id)
//     {
//         ref_id = r_id;
//     }
// 
//     void set_start_loc(size_t s_loc)
//     {
//         start_loc = s_loc;
//     }
// 
//     void set_end_loc(size_t e_loc)
//     {
//         end_loc = e_loc;
//     }
// };
// 
// vector<SegInfo> seginfo_vec;
// 
// class Consensus: public loon::Consensus
// {
// public:
//     vector<size_t> path;
// public:
//     Consensus()
//     {
//     }
// 
//     Consensus(const loon::Consensus& base): loon::Consensus(base)
//     {
//     }
// 
//     void append_path(size_t id) // the vertex id of the breakpoint graph
//     {
//         path.push_back( id );
//     }
// 
//     void print_path(ostream& out, size_t id)
//     {
//         out << id << ' ' << path.size() << endl;
//         for(vector<size_t>::iterator p_it = path.begin(); p_it != path.end(); ++p_it)
//         {
//             size_t seg_id = *p_it;
//             bool is_backward = seg_id & 1;
//             seg_id >>= 1;
//             out << seginfo_vec[ seg_id ].ref_id << ' '
//                 << seginfo_vec[ seg_id ].start_loc << ' '
//                 << seginfo_vec[ seg_id ].end_loc << ' '
//                 << is_backward << endl;
//         }
//     }
// };
// 
// template<class HeapArrayType>
// class MaxHeap
// {
// protected:
//     const HeapArrayType& a;
//     vector<size_t> heap;
// public:
//     MaxHeap(const HeapArrayType& aa): a(aa)
//     {
//         size_t n = a.size();
//         if(n == 0)
//             return;
//         heap.resize(n);
//         // init
//         for(size_t i = 0; i < n; ++i)
//             heap[i] = i;
//         //make heap
//         for(size_t i = (n-1)>> 1; i != -1; --i)
//             moveDown(i);
//     }
// 
//     void moveDown(size_t parent)
//     {
//         size_t n = heap.size();
//         if(n <= 1)  return;
//         size_t max_child = parent << 1;
//         while(max_child < n)
//         {
//             if((max_child | 1) < n && a[ heap[max_child | 1] ] > a[ heap[max_child] ])
//                 max_child |= 1;
//             if(a[ heap[max_child >> 1] ] < a[ heap[max_child] ])
//                 swap(heap[max_child >> 1], heap[max_child]);
//             else    return;
//             max_child <<= 1;
//         }
//     }
// 
//     size_t getMaxId() const
//     {
//         return heap.empty() ? -1 : heap[0];
//     }
// 
//     void removeMax()
//     {
//         swap(heap[0], heap.back());
//         heap.pop_back();
//         moveDown(0);
//     }
// };

class ArcPair
{
public:
    size_t source, target;
public:
    ArcPair(size_t s = 0, size_t t = 0): source(s), target(t)
    {
    }
    bool operator<(const ArcPair& t) const
    {
        return (source < t.source || (source == t.source && target < t.target));
    }
};

class PrimalGraph
{
public:
    vector<size_t> vid_to_newid;
    vector<size_t> newid_to_vid;
    vector<vector<size_t> > edge_list;
    vector<vector<long long> > source_cutoffs, target_cutoffs;
    //vector<ArcPair> edge_endpoints;
    map<ArcPair, size_t> edge_set;

    vector<long long> edge_weight;
    vector<size_t> best_edge_id;
    vector<long long> best_source_cutoffs, best_target_cutoffs;
public:
    PrimalGraph(size_t n = 0): 
        vid_to_newid(n << 1, loon::INF)
    {}

    bool is_unmapped_vertex(size_t id)
    {
        return (vid_to_newid[id] == loon::INF);
    }

    size_t get_node_size() const
    {
        return newid_to_vid.size();
    }

    size_t get_edge_size() const
    {
        return edge_set.size();
    }

    size_t get_original_vid(size_t id) const
    {
        return newid_to_vid[ id ];
    }

    size_t get_new_vid(size_t id) const
    {
        return vid_to_newid[ id ];
    }

    bool get_edge_string_info(size_t ori_u, size_t ori_v, 
            size_t& e_id, long long& left_cutoff, long long& right_cutoff)
    {
        bool forward_edge = false;
        if(ori_u < ori_v)
        {
            ArcPair arc(ori_u, ori_v);
            size_t edge_id = edge_set[ arc ];
            e_id = best_edge_id[ edge_id ];
            left_cutoff = best_source_cutoffs[ edge_id ];
            right_cutoff = best_target_cutoffs[ edge_id ];

            //if(ori_u & 1)// u_t: (u_t, v_t, e_r) or (u_t, v_h, e_r)
            //    forward_edge = e_id & 1;// e_id & 1 ? true : false
            //else// u_h: (u_h, v_t, e_f) or (u_h, v_h, e_f)
            //    forward_edge = !(e_id & 1); // e_id & 1 ? false : true

            forward_edge = !((ori_u ^ e_id) & 1);
        }
        else
        {
            ArcPair arc(ori_v, ori_u);
            size_t edge_id = edge_set[ arc ];
            e_id = best_edge_id[ edge_id ];
            left_cutoff = best_target_cutoffs[ edge_id ];
            right_cutoff = best_source_cutoffs[ edge_id ];

            //if(ori_v & 1)// v_t: (v_t, u_h, e_f) or (v_t, u_t, e_f)
            //    forward_edge = !(e_id & 1); // e_id & 1 ? false : true;
            //else// v_h: (v_h, u_h, e_r) or (v_h, u_t, e_r)
            //    forward_edge = e_id & 1; // e_id & 1 ? true : false;
            forward_edge = (ori_v ^ e_id) & 1;
        }
        e_id >>= 1;
        return forward_edge;
    }

    void add_new_vertex(size_t v)
    {
        vid_to_newid[ v ] = newid_to_vid.size();
        newid_to_vid.push_back( v );
    }

    void add_edge(size_t source, size_t target, long long source_cutoff, long long target_cutoff, size_t e_id)
    {
        ArcPair arc(source, target);
        size_t id;
        map<ArcPair, size_t>::iterator it = edge_set.find( arc );
        if( it == edge_set.end())
        {
            id = edge_list.size();
            edge_set[ arc ] = id;
            //edge_endpoints.push_back( arc );

            edge_list.push_back( vector<size_t>() );
            source_cutoffs.push_back( vector<long long>() );
            target_cutoffs.push_back( vector<long long>() );
        }
        else
            id = it->second;
        edge_list[ id ].push_back( e_id );
        source_cutoffs[ id ].push_back( source_cutoff );
        target_cutoffs[ id ].push_back( target_cutoff );

        if( get_new_vid(source) == loon::INF )
        {
            add_new_vertex(source);
            add_new_vertex(source ^ 1);
        }
        if( get_new_vid(target) == loon::INF )
        {
            add_new_vertex(target);
            add_new_vertex(target ^ 1);
        }
    }

    // u, v, e_id are indices of the original vertices, and edges
    void add_edge_fl_fr(size_t u, size_t v, long long uside_cutoff, long long vside_cutoff, size_t e_id)// case 1
    {
        if(u == v)  return;
        if(u < v)
            add_edge(head_idx(u), tail_idx(v), uside_cutoff, vside_cutoff, forward_eid(e_id));
        else
            add_edge(tail_idx(v), head_idx(u), vside_cutoff, uside_cutoff, forward_eid(e_id));
    }

    void add_edge_fl_rl(size_t u, size_t v, long long uside_cutoff, long long vside_cutoff, size_t e_id)// case 2
    {
        if(u == v)  return;
        if(u < v)
            add_edge(head_idx(u), head_idx(v), uside_cutoff, vside_cutoff, forward_eid(e_id));
        else
            add_edge(head_idx(v), head_idx(u), vside_cutoff, uside_cutoff, reverse_eid(e_id));
    }

    void add_edge_rr_fr(size_t u, size_t v, long long uside_cutoff, long long vside_cutoff, size_t e_id)// case 3
    {
        if(u == v)  return;
        if(u < v)
            add_edge(tail_idx(u), tail_idx(v), uside_cutoff, vside_cutoff, reverse_eid(e_id));
        else
            add_edge(tail_idx(v), tail_idx(u), vside_cutoff, uside_cutoff, forward_eid(e_id));
    }

    void add_edge_rr_rl(size_t u, size_t v, long long uside_cutoff, long long vside_cutoff, size_t e_id)// case 4
    {
        if(u == v)  return;
        if(u < v)
            add_edge(tail_idx(u), head_idx(v), uside_cutoff, vside_cutoff, reverse_eid(e_id));
        else
            add_edge(head_idx(v), tail_idx(u), vside_cutoff, uside_cutoff, reverse_eid(e_id));
    }

    void add_edges(size_t e_id, 
            const vector<size_t>& u_vertices, const vector<size_t>& v_vertices, 
            const vector<long long>& left_cutoffs, const vector<long long>& right_cutoffs,
            int edge_case)
    {
        if(u_vertices.empty() || v_vertices.empty())
            return;
        size_t u_n = u_vertices.size();
        size_t v_n = v_vertices.size();

        for(size_t u_i = 0; u_i < u_n; ++u_i)
        {
            for(size_t v_i = 0; v_i < v_n; ++v_i)
            {
                if(edge_case == 1)
                    add_edge_fl_fr(u_vertices[u_i], v_vertices[v_i], left_cutoffs[u_i], right_cutoffs[v_i], e_id);
                else if(edge_case == 2)
                    add_edge_fl_rl(u_vertices[u_i], v_vertices[v_i], left_cutoffs[u_i], right_cutoffs[v_i], e_id);
                else if(edge_case == 3)
                    add_edge_rr_fr(u_vertices[u_i], v_vertices[v_i], left_cutoffs[u_i], right_cutoffs[v_i], e_id);
                else
                    add_edge_rr_rl(u_vertices[u_i], v_vertices[v_i], left_cutoffs[u_i], right_cutoffs[v_i], e_id);
            }
        }
    }

    void add_edges(loon::PathCover& path_cover)
    {
        vector<size_t> deg(vid_to_newid.size(), 0);
        for(map<ArcPair, size_t>::iterator it = edge_set.begin();
                it != edge_set.end(); ++it)
        {
            ++deg[ get_new_vid( it->first.source ) ];
            ++deg[ get_new_vid( it->first.target ) ];
        }
        for(map<ArcPair, size_t>::iterator it = edge_set.begin();
                it != edge_set.end(); ++it)
        {
            size_t new_u = get_new_vid( it->first.source);
            size_t new_v = get_new_vid( it->first.target);
            double multi_deg = double(deg[ new_u ]) * deg[ new_v ];
            path_cover.addEdge( new_u,
                                new_v,
                                edge_weight[ it->second ] / multi_deg);

        }
        /// for(map<ArcPair, size_t>::iterator it = edge_set.begin();
        ///         it != edge_set.end(); ++it)
        /// {
        ///     //if(edge_weight[ it->second ] >= 15)

        ///     path_cover.addEdge( get_new_vid( it->first.source),
        ///                         get_new_vid( it->first.target),
        ///                         edge_weight[ it->second ]);

        ///     // cerr <<  get_new_vid( it->first.source) << ' '
        ///     //      << get_new_vid( it->first.target) << ' '
        ///     //      << edge_weight[ it->second ] << endl;
        /// }
    }

    bool check_connectivity(string& sleft, string& sright, string& e, loon::Aligner& verify_aligner)
    {
        verify_aligner.set_strs(sleft, e);
        if(!( (run_options & 0x10) ? 
                    verify_aligner.align_right(loon::blasr_align_choice, true) :
                    verify_aligner.align_right(true)))
            return false;
        verify_aligner.set_strs(e, sright);
        if(!( (run_options & 0x10) ? 
                    verify_aligner.align_right(loon::blasr_align_choice, true) :
                    verify_aligner.align_right(true)))
            return false;
        return true;
    }

    void merge_similar_edges(loon::Multiseq& long_reads)
    {
        const size_t nn = edge_list.size();
        edge_weight.assign(nn, 0);
        best_edge_id.assign(nn, loon::INF);
        best_source_cutoffs.assign(nn, -1);
        best_target_cutoffs.assign(nn, -1);

        for(size_t i = 0; i < nn; ++i)
        {
            size_t i_size = edge_list[i].size();
            if(i_size == 1)
            {
                best_edge_id[ i ] = 0;
                edge_weight[ i ] = 1;
                continue;
            }

            string infile = loon::mhap_get_infile_name( i );
            ofstream fout( infile.c_str() );
            for(size_t j = 0; j < i_size; ++j)
            {
                fout << ">" << j << endl;
                if(edge_list[i][j] & 1)
                    fout << long_reads[ edge_list[i][j] >> 1 ].rev << endl;
                else
                    fout << long_reads[ edge_list[i][j] >> 1 ].seq << endl;
            }
            fout.close();
        }
        // run mhap in parallel
        loon::mhap_pick_best_edge_id(nn, best_edge_id, edge_weight);

        for(size_t i = 0; i < nn; ++i)
            if(best_edge_id[ i ] != loon::INF)
            {
                best_source_cutoffs[ i ] = source_cutoffs[i][ best_edge_id[i] ];
                best_target_cutoffs[ i ] = target_cutoffs[i][ best_edge_id[i] ];
                best_edge_id[i] = edge_list[ i ][ best_edge_id[i] ];
            }
        edge_list.clear();
        source_cutoffs.clear();
        target_cutoffs.clear();
    }
};

bool add_edge(vector<size_t>& left_vertices, vector<size_t>& right_vertices,
        vector<long long>& left_boundary, vector<long long>& right_boundary,
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
    {
        right_vertices.push_back( str2lu( delta.ref_tag ) );
        right_boundary.push_back( delta.qry_len - qry_left_ext + ref_left_ext );
    }
    else
    {
        left_vertices.push_back( str2lu(delta.ref_tag) );
        left_boundary.push_back( delta.clusters.back().qry_end + ref_right_ext );
    }
    return true;
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

void do_consensus_save_edge(PrimalGraph& p_graph, size_t u_id, size_t v_id,
        vector<string>& contigs, vector<string>& contigs_score,
        loon::Multiseq& tmp_ref, vector<string>& ctg_score_tmp,
        loon::Multiseq& long_reads)
{
    size_t e_id;
    long long left_cutoff, right_cutoff;
    bool forward_edge = p_graph.get_edge_string_info( u_id, v_id, e_id, left_cutoff, right_cutoff );
    if(left_cutoff + right_cutoff > long_reads[ e_id ].seq.length()) // u and v overlap
    {
        size_t overlap_len = left_cutoff + right_cutoff - long_reads[ e_id ].seq.length();
        if(v_id & 1) // forward 
        {
            contigs.back().append( tmp_ref[ v_id>>1 ].seq.begin() + overlap_len, tmp_ref[ v_id >>1 ].seq.end() );
            contigs_score.back().append( ctg_score_tmp[ v_id>>1 ].begin() + overlap_len, ctg_score_tmp[ v_id>>1 ].end() );
        }
        else
        {
            string tmp_str(tmp_ref[ v_id>>1 ].seq.begin() , tmp_ref[ v_id >>1 ].seq.end() - overlap_len);
            loon::compute_reverse_complement( tmp_str );
            contigs.back().append( tmp_str );
            contigs_score.back().append( ctg_score_tmp[ v_id >> 1 ].rbegin() + overlap_len, ctg_score_tmp[ v_id>>1 ].rend() );
        }
    }
    else
    {// u and v have gap
        // fill the gap
        size_t gap_len = long_reads[ e_id ].seq.length() - left_cutoff - right_cutoff;
        if(forward_edge)
            contigs.back().append( long_reads[ e_id ].seq.begin() + left_cutoff, long_reads[ e_id ].seq.end() - right_cutoff );
        else
            contigs.back().append( long_reads[ e_id ].rev.begin() + left_cutoff, long_reads[ e_id ].rev.end() - right_cutoff );
        contigs_score.back().append( gap_len, '!' );

        // append v
        if(v_id & 1) // forward
        {
            contigs.back().append( tmp_ref[ v_id>>1 ].seq );
            contigs_score.back().append( ctg_score_tmp[ v_id>>1 ]);
        }
        else
        {
            loon::compute_reverse_complement( tmp_ref[ v_id >> 1 ].seq );
            contigs.back().append( tmp_ref[v_id >> 1].seq );
            contigs_score.back().append( ctg_score_tmp[ v_id >> 1 ].rbegin(), ctg_score_tmp[ v_id>>1 ].rend() );
        }
    }
    tmp_ref[ v_id >> 1].seq.clear();
    ctg_score_tmp[ v_id >> 1].clear();
}

void paths_refinement(vector<string>& contigs, vector<string>& contigs_score, loon::Multiseq& long_reads, vector<bool>& used_long_reads)
{
    //keep_chosen_vector(long_reads, used_long_reads);
    used_long_reads.clear();
    size_t contigs_size = contigs.size();
    loon::Multiseq tmp_ref;
    vector<string> ctg_score_tmp;
    tmp_ref.resize( contigs_size );

    for(size_t i = 0; i < contigs_size; ++i)
        tmp_ref[ i ].seq.swap( contigs[ i ] );
    ctg_score_tmp.swap( contigs_score );
    contigs.clear();    contigs.reserve( contigs_size );
    contigs_score.clear();  contigs_score.reserve( contigs.size() );

    nucmer_Deltas deltas_forward, deltas_backward;
    loon::all_v_all_matching(tmp_ref, long_reads, deltas_forward, deltas_backward);
    sort(deltas_forward.begin(), deltas_forward.end(), deltas_tag_cmp);
    sort(deltas_backward.begin(), deltas_backward.end(), deltas_tag_cmp);
    size_t fi = 0, bi = 0;
    size_t fn = deltas_forward.size(), bn = deltas_backward.size();
    string ref_tag;
    INFO_PRINT("start path quality refinement");
    while(fi < fn || bi < bn)
    {
        loon::Aligner alner;
        loon::MSA_Consensus msa_consensus;
        // update ref_tag
        if(fi < fn)
        { 
            ref_tag = deltas_forward[ fi ].ref_tag;
            if(bi < bn && ref_tag > deltas_backward[ bi ].ref_tag)
                ref_tag = deltas_backward[ bi ].ref_tag;
        }
        else    ref_tag = deltas_backward[ bi ].ref_tag;
        size_t ref_id = str2lu( ref_tag );
        msa_consensus.add_ref( &tmp_ref[ ref_id ].seq, &ctg_score_tmp[ ref_id ] );
        alner.set_ref( tmp_ref[ ref_id ].seq );
        for(; fi < fn && deltas_forward[ fi ].ref_tag == ref_tag; ++fi)
        {
            size_t qry_id = str2lu( deltas_forward[fi].qry_tag );
            alner.set_qry( long_reads[ qry_id ].seq );
            if(alner.fill_gaps(false, true, deltas_forward[ fi ]))
            {
                msa_consensus.add_alignment( alner.XX_start - 1, alner.XX_end,
                        alner.YY_start - 1, alner.YY_end,
                        long_reads[ qry_id ].seq, alner.indels);
            }
        }
        for(; bi < bn && deltas_backward[ bi ].ref_tag == ref_tag; ++bi )
        {
            size_t qry_id = str2lu( deltas_backward[bi].qry_tag );
            alner.set_qry( long_reads[ qry_id ].rev );
            if(alner.fill_gaps(false, true, deltas_backward[ bi ]))
            {
                msa_consensus.add_alignment( alner.XX_start - 1, alner.XX_end,
                        alner.YY_start - 1, alner.YY_end,
                        long_reads[ qry_id ].rev, alner.indels);
            }
        }
        msa_consensus.sort();
        msa_consensus.do_consensus( true );
        vector<string>& msa_cons = msa_consensus.get_consensus();
        vector<string>& msa_qual = msa_consensus.get_consensus_score();
        size_t msa_n = msa_cons.size();
        for(size_t msa_i = 0; msa_i < msa_n; ++msa_i)
        {
            contigs.push_back( string() );
            contigs.back().swap( msa_cons[ msa_i ] );
            contigs_score.push_back( string() );
            contigs_score.back().swap( msa_qual[ msa_i ] );
        }
    }
}

void do_bridging(vector<string>& contigs, vector<string>& contigs_score, loon::Multiseq& long_reads, vector<size_t>& lr_idx)
{
    if(contigs.size() <= 1)
        return;
    if(consensus_strategy == 2)
        keep_chosen_vector(long_reads, lr_idx);
    lr_idx.clear();

    vector<bool> used_long_reads;

    const size_t lr_n = long_reads.size();
    used_long_reads.assign(lr_n, false);

    size_t contigs_size = contigs.size();
    INFO_PRINT("Convert %lu contigs to strings", contigs_size);
    vector<string> ctg_score_tmp;
    loon::Multiseq tmp_ref;
    tmp_ref.resize( contigs_size );
    ctg_score_tmp.swap( contigs_score );
    for(size_t j = 0; j < contigs_size; ++j)
        tmp_ref[j].seq.swap( contigs[j] ); // use swap to save memory and running time

    INFO_PRINT("All contigs vs. all long reads matching");
    nucmer_Deltas deltas_forward, deltas_backward;
    loon::all_v_all_matching(tmp_ref, long_reads, 
            deltas_forward, deltas_backward, true);

    vector<vector<size_t> > fleft_vertices( lr_n ), fright_vertices( lr_n ),
                            rleft_vertices( lr_n ), rright_vertices( lr_n );
    vector<vector<long long> > fleft_boundary( lr_n ), fright_boundary( lr_n ),
                              rleft_boundary( lr_n ), rright_boundary( lr_n );

    INFO_PRINT("Compute vertices for edges");
    progress.start( deltas_forward.size() + deltas_backward.size() );
    for(size_t ii = 0; ii < deltas_forward.size(); ++ii)
    {
        size_t lr_i = str2lu(deltas_forward[ ii ].qry_tag);
        used_long_reads[ lr_i ] = add_edge(fleft_vertices[ lr_i ], fright_vertices[ lr_i ],
                                           fleft_boundary[ lr_i ], fright_boundary[ lr_i ],
                                            deltas_forward[ ii ]) ||
                                    used_long_reads[ lr_i ];
        progress.progress(1);
    }
    for(size_t ii = 0; ii < deltas_backward.size(); ++ii)
    {
        size_t lr_i = str2lu(deltas_backward[ ii ].qry_tag);
        used_long_reads[ lr_i ] = add_edge(rleft_vertices[ lr_i ], rright_vertices[ lr_i ],
                                            rleft_boundary[ lr_i ], rright_boundary[ lr_i ], 
                                            deltas_backward[ ii ]) ||
                                    used_long_reads[ lr_i ];
        progress.progress(1);
    }
    progress.stop();
    
    INFO_PRINT("Build the primal graph");
    PrimalGraph p_graph( contigs_size );
    for(size_t lr_i = 0; lr_i < lr_n; ++lr_i)
    {
        // case 1: (fl, fr)
        p_graph.add_edges(lr_i, fleft_vertices[lr_i], fright_vertices[lr_i], 
                fleft_boundary[lr_i], fright_boundary[lr_i], 1);
        // case 2: (fl, rl)
        p_graph.add_edges(lr_i, fleft_vertices[lr_i], rleft_vertices[lr_i], 
                fleft_boundary[lr_i], rleft_boundary[lr_i], 2);
        // case 3: (rr, fr)
        p_graph.add_edges(lr_i, rright_vertices[lr_i], fright_vertices[lr_i], 
                rright_boundary[lr_i], fright_boundary[lr_i], 3);
        // case 4: (rr, rl)
        p_graph.add_edges(lr_i, rright_vertices[lr_i], rleft_vertices[lr_i],
                rright_boundary[lr_i], rleft_boundary[lr_i], 4);

    }
    fleft_vertices.clear();     fleft_boundary.clear();
    fright_vertices.clear();    fright_boundary.clear();
    rleft_vertices.clear();     rleft_boundary.clear();
    rright_vertices.clear();    rright_boundary.clear();

    INFO_PRINT("Merge similar edges");
    p_graph.merge_similar_edges(long_reads);

    INFO_PRINT("Build actual graph");
    loon::PathCover path_cover;
    size_t m = p_graph.get_edge_size();
    path_cover.reserveNode( p_graph.get_node_size() );
    path_cover.reserveEdge( m );
    path_cover.genNodes();

    //DERR << "Print the graph" << endl;
    //cerr << p_graph.get_node_size() << ' ' << m << endl;
    p_graph.add_edges( path_cover );

    INFO_PRINT("run path cover");
    path_cover.run();

    INFO_PRINT("Consensus paths: [%lu paths]", path_cover.paths.size());
    contigs.clear();
    for(size_t i = 0; i < path_cover.paths.size(); ++i)
    {
        if(path_cover.paths[i].path_type == 1)
        {
            size_t u_id = p_graph.get_original_vid( path_cover.paths[i].edges[0] ) >> 1;
            contigs.push_back( string() );
            contigs.back().swap( tmp_ref[ u_id ].seq );
            contigs_score.push_back( string() );
            contigs_score.back().swap( ctg_score_tmp[ u_id ] );

        }
        else
        {
            size_t u_id, v_id;

            u_id = p_graph.get_original_vid( path_cover.paths[i].edges[0] );
            contigs.push_back( string() );
            contigs_score.push_back( string() );
            contigs.back().swap( tmp_ref[ u_id >> 1 ].seq );
            contigs_score.back().swap( ctg_score_tmp[ u_id >>1 ] );

            if(u_id & 1)
            {
                loon::compute_reverse_complement( contigs.back() );
                reverse( contigs_score.back().begin(), contigs_score.back().end());
            }

            for(size_t j = 0; j < path_cover.paths[i].edges.size(); j+=2)
            {
                u_id = p_graph.get_original_vid( path_cover.paths[i].edges[j] );
                v_id = p_graph.get_original_vid( path_cover.paths[i].edges[j+1] );

                do_consensus_save_edge(p_graph, u_id, v_id, contigs, contigs_score,
                        tmp_ref, ctg_score_tmp, long_reads);
            }
        }
    }
    deltas_forward.clear(); deltas_backward.clear();

    // path quality refinement
    INFO_PRINT("last alignment for path quality refinement");
    paths_refinement(contigs, contigs_score, long_reads, used_long_reads);
    // add unaligned contigs
    size_t lol_cnt = 0;
    for(size_t i = 0; i < tmp_ref.size(); ++i)
    {
        if(!tmp_ref[ i ].seq.empty())
        {
            ++lol_cnt;
            contigs.push_back( string() );
            contigs_score.push_back( string() );
            contigs.back().swap( tmp_ref[ i ].seq );
            contigs_score.back().swap( ctg_score_tmp[ i ] );
        }
    }
    INFO_PRINT("# of unaligned contigs: %lu", lol_cnt);
}

/*================== class NCAlignment ==================
========================= start =======================*/
class NCAlignment
{
public:
    size_t ref_start, ref_end;
    string* strp;
    tigrinc::Nucmer_Delta* delta;
public:
    NCAlignment(size_t rs = 0, size_t rt = 0, string* sp = NULL, tigrinc::Nucmer_Delta* d = NULL):
        ref_start(rs), ref_end(rt), strp(sp), delta(d)
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
    void add_coverage(tigrinc::Nucmer_Delta& d, string* strp);
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

void Regions::add_coverage(tigrinc::Nucmer_Delta& d, string* strp)
{
    if(d.clusters.empty())
        return;

    Region res;
    res.start = d.clusters.front().ref_start;
    res.end = d.clusters.back().ref_end;
    res.long_reads = new list<NCAlignment>();
    res.long_reads->push_back( NCAlignment(res.start, res.end, strp, &d) );

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


void read_contigs(vector<string>& contigs, vector<string>& contigs_score)
{
    ifstream fin(ctg_file.c_str());
    loon::check_file_open( fin, ctg_file );

    string line;
    while(fin >> line)
    {
        contigs.push_back( string() );
        contigs.back().swap( line );
        contigs_score.push_back( string() );
        fin >> contigs_score.back();
    }

    fin.close();

    //string segfile = prefix + ".dis";
    //fin.open(segfile.c_str());
    //loon::check_file_open(fin, segfile);

    //SegInfo tmp_seginfo;
    //while(fin >> tmp_seginfo.ref_id >> tmp_seginfo.start_loc >> tmp_seginfo.end_loc)
    //    seginfo_vec.push_back(tmp_seginfo);

    //if(seginfo_vec.size() != cnt)
    //{
    //    ERROR_PRINT("[%s] and [%s] don\'t match", segfile.c_str(), ctg_file.c_str());
    //    exit(1);
    //}

    //fin.close();
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

void save_final_result(vector<string>& contigs, vector<string>& contigs_score)
{
    ofstream fout_finalseq, fout_finalfasta, fout_path;
    open_file( ".final.seq", fout_finalseq );
    open_file( ".final.fasta", fout_finalfasta);

    loon::Fasta tmp_fasta;
    size_t contig_len = 0;
    for(size_t i = 0; i < contigs.size(); ++i)
    {
        fout_finalseq << contigs[i] << endl;
        fout_finalseq << contigs_score[i] << endl;

        contig_len += contigs[i].length();

        ostringstream oss;
        oss << "cons_" << i;
        tmp_fasta.tag = oss.str();
        tmp_fasta.seq.swap( contigs[i] );
        tmp_fasta.save_fasta( fout_finalfasta );
    }
    
    fout_finalseq.close();
    fout_finalfasta.close();

    INFO_PRINT("%lu contigs, total length = %lu", contigs.size(), contig_len);

#ifdef DEBUG_CHGCNT
    INFO_PRINT("%lu changes", chgcnt - contigs.size());
#endif
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
            regions.add_coverage( *df, &qry[ qry_id ].seq);
            qry_used[ qry_id ] = true;
            ref_used[ str2lu(df->ref_tag) ] = true;
        }
    }
    else
    {
        if(b_score >= min_score)
        {
            size_t qry_id = str2lu( db->qry_tag );
            regions.add_coverage(*db, &qry[ qry_id ].rev);
            qry_used[ qry_id ] = true;
            ref_used[ str2lu(db->ref_tag) ] = true;
        }
    }
}

void consensus_regions(Regions& regions, vector<string>& contigs, vector<string>& contigs_score,
        ostream& fout_cons_fasta, ostream& fout_cons_seq, string& ref_seq, size_t ref_id)
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
            //seg_id = seginfo_vec.size();
            start_loc = rgn_it->start;
            end_loc = rgn_it->end;
            //seginfo_vec.push_back( SegInfo(ref_id, start_loc, end_loc) );

            contigs.push_back( ref_seq.substr(start_loc - 1, end_loc - start_loc + 1) );
            contigs_score.push_back( string(end_loc - start_loc + 1, '!') );
        }
        else
        {
            // new efficient implementation for consensus
            loon::Aligner alner;
            loon::MSA_Consensus msa_consensus;
            msa_consensus.add_ref(&ref_seq, NULL);
            alner.set_ref(ref_seq);
            for(list<NCAlignment>::iterator nc_it = rgn_it -> long_reads -> begin();
                    nc_it != rgn_it -> long_reads ->end(); ++nc_it)
            {
                alner.set_qry( *(nc_it->strp) );
                if(alner.fill_gaps(false, true, *(nc_it->delta)))
                {
                    msa_consensus.add_alignment( alner.XX_start - 1, alner.XX_end, alner.YY_start - 1, alner.YY_end,
                            *(nc_it->strp), alner.indels);
                }
            }
            msa_consensus.sort();
            msa_consensus.do_consensus();
            vector<string>& msa_cons = msa_consensus.get_consensus();
            vector<string>& msa_qual = msa_consensus.get_consensus_score();
            size_t msa_n = msa_cons.size();
            for(size_t msa_i = 0; msa_i < msa_n; ++msa_i)
            { 
                if(run_options & 1)
                {
                    fout_cons_seq << msa_cons[ msa_i ] << endl;
                    fout_cons_seq << msa_qual[ msa_i ] << endl;

                    tmp_fasta.seq = msa_cons[ msa_i ];
                    ostringstream oss;
                    oss << "cons_" << contigs.size();
                    tmp_fasta.tag = oss.str();
                    tmp_fasta.save_fasta( fout_cons_fasta );
                }
                contigs.push_back( string() );          contigs.back().swap( msa_cons[ msa_i ] );
                contigs_score.push_back( string() );    contigs_score.back().swap( msa_qual[ msa_i ] );
            }
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
        vector<string>& contigs, vector<string>& contigs_score, vector<size_t>& lr_idx)
{
    size_t should_have_regions_cnt = 0;
    ofstream fout_cons_fasta, fout_cons_seq, fout_idx;
    if(run_options & 1)
    {
        open_file(".cons.fasta", fout_cons_fasta);
        open_file(".cons.seq", fout_cons_seq);
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
            consensus_regions(regions, contigs, contigs_score, fout_cons_fasta, fout_cons_seq, ref[ ref_id ].seq, ref_id);
            last_tag = min(deltas_forward[i].ref_tag, deltas_backward[j].ref_tag);
        }
    }

    while(i < deltas_fn)
    {
        if(deltas_forward[i].ref_tag != last_tag)
        {
            should_have_regions_cnt += regions.size();
            size_t ref_id = str2lu( last_tag );
            consensus_regions(regions, contigs, contigs_score, fout_cons_fasta, fout_cons_seq, ref[ ref_id ].seq, ref_id);
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
            consensus_regions(regions, contigs, contigs_score, fout_cons_fasta, fout_cons_seq, ref[ ref_id ].seq, ref_id);
            last_tag = deltas_backward[j].ref_tag;
        }
        try_add_regions(NULL, &deltas_backward[j], qry, regions, ref_used, qry_used);
        ++j;
        progress.progress(1);
    }
    progress.stop();
    should_have_regions_cnt += regions.size();
    size_t ref_id = str2lu(last_tag);
    consensus_regions(regions, contigs, contigs_score, fout_cons_fasta, fout_cons_seq, ref[ ref_id ].seq, ref_id);
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
                contigs.push_back( ref[i].seq );
                contigs_score.push_back( string(ref[i].seq.length(), '!') );

                //seginfo_vec.push_back( SegInfo(i, 1, ref[i].seq.length()) );

                if(run_options & 1)
                {
                    fout_cons_seq << contigs.back() << endl;
                    fout_cons_seq << contigs_score.back() << endl;
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
        fout_cons_seq.close();
        fout_idx.close();
    }

    // save disassembling info
    //ofstream fout_seginfo;
    //open_file(".dis", fout_seginfo);
    //for(size_t i = 0; i < seginfo_vec.size(); ++i)
    //{
    //    fout_seginfo << seginfo_vec[i].ref_id << ' ' 
    //                << seginfo_vec[i].start_loc << ' ' 
    //                << seginfo_vec[i].end_loc << endl;
    //}
    //fout_seginfo.close();
}

void help_msg(char* exec)
{
    cerr << "Version: 2.2.1 !!!" << endl;
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
    cerr << "\t                     * 1: do consensus (faster)" << endl;
    cerr << "\t                     * 2: do consensus (higher quality)." << endl;
    cerr << "\t-B               Don't bridge the unaligned references" << endl;
    cerr << "\t-j threads       # of threads" << endl;
    cerr << "\t-h               show help info" << endl;
    cerr << "Notice:" << endl;
    cerr << "\tThe -c and -i works only when only do bridging" << endl;
    cerr << "\tFor lower coverage data, you probably want to use -C 2 to get higher quality results" << endl;
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

    if(((run_options & 3) == 2) && (ctg_file.empty() || idx_file.empty()))
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

    vector<string> contigs, contigs_score;
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
        do_alignment_and_consensus(ref, qry, contigs, contigs_score, lr_idx);
    }
    
    if(run_options & 6)
    { // do bridging step
        if((run_options & 3) == 2)
        { // read contigs file and index file
            INFO_PRINT("read contigs file and index file");

            read_contigs(contigs, contigs_score);
            read_idx(lr_idx);
        }
        
        INFO_PRINT("do bridging");
        do_bridging(contigs, contigs_score, qry, lr_idx);

        INFO_PRINT("save bridging result");
        save_final_result(contigs, contigs_score);
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
