#ifndef __ALIGNER_H
#define __ALIGNER_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib>

#include "baseheader.h"
#include "nucmer.h"
#include "util.h"
#include "multi-array.h"


namespace loon
{

const size_t max_ext = 50;
//const size_t default_force_max_diff = 200;
const size_t default_force_max_diff = 100;
const size_t max_end_padding = 1000;

class Aligner
{
public:
    const std::string* XX;
    const std::string* YY;

    int error_code;
    bool simple_align;

    size_t max_gap;
    size_t max_diff;

    Array2D<double> a;

    size_t XX_start, XX_end;
    size_t YY_start, YY_end;
    std::vector<long> indels;

    const void* cmp_obj;
    double (*mismatch)(const void*, size_t, size_t);
    double (*delta)(bool);
    size_t uabsminus(size_t a, size_t b);
public:
    Aligner(size_t max_size = 4096);
    void set_simple_align();
    void unset_simple_align();
    void set_strs(const string& a, const string& b);
    void set_ref(const string& a);
    void set_qry(const string& b);
    void set_cmp_obj();
    void set_cmp_obj(const void* obj);
    void set_parameters(size_t max_gap = 20, size_t max_diff = 10);

    /* if force_align is true, max_gap will be unused, and the default 
       max_diff will be much larger, unless a even larger max_diff is 
       set by calling set_parameters() */
    bool align_full(void (*compute_deltas_func)(const std::string&, 
                    const std::string&, nucmer_Deltas&, size_t, int),
            bool force_align = false,
            double (*mismatch)(const void*, size_t, size_t) = default_mismatch,
            double (*delta)(bool) = default_delta);
    bool align_full(bool force_align = false,
            double (*mismatch)(const void*, size_t, size_t) = default_mismatch,
            double (*delta)(bool) = default_delta);
    bool align_full(Suffixtree& stree,
            Uchar* ref,
            Uint ref_len,
            bool force_align = false,
            double (*mismatch)(const void*, size_t, size_t) = default_mismatch,
            double (*delta)(bool) = default_delta);

    // for align_left and align_right, the code for the case that the suffix tree 
    // is provided is not correct, so don't call those functions
    bool align_left(void (*compute_deltas_func)(const std::string&,
                    const std::string&, nucmer_Deltas&, size_t, int),
            bool force_align = false,
            size_t rightmost_pos = INF,
            double (*mismatch)(const void*, size_t, size_t) = default_mismatch,
            double (*delta)(bool) = default_delta);
    bool align_left(bool force_align = false,
            size_t rightmost_pos = INF,
            double (*mismatch)(const void*, size_t, size_t) = default_mismatch,
            double (*delta)(bool) = default_delta);
    bool align_left(Suffixtree& stree,
            Uchar* ref,
            Uint ref_len,
            bool force_align = false,
            size_t rightmost_pos = INF,
            double (*mismatch)(const void*, size_t, size_t) = default_mismatch,
            double (*delta)(bool) = default_delta);

    bool align_right(void (*compute_deltas_func)(const std::string&,
                    const std::string&, nucmer_Deltas&, size_t, int),
            bool force_align = false,
            size_t leftmost_pos = INF,
            double (*mismatch)(const void*, size_t, size_t) = default_mismatch,
            double (*delta)(bool) = default_delta);
    bool align_right(bool force_align = false,
            size_t leftmost_pos = INF,
            double (*mismatch)(const void*, size_t, size_t) = default_mismatch,
            double (*delta)(bool) = default_delta);
    bool align_right(Suffixtree& stree,
            Uchar* ref,
            Uint ref_len,
            bool force_align = false,
            size_t leftmost_pos = INF,
            double (*mismatch)(const void*, size_t, size_t) = default_mismatch,
            double (*delta)(bool) = default_delta);

    /*---------------- This is the dividing line --------------------*/

    bool simple_align_full_primal(bool force_align,
            Suffixtree* stree,
            Uchar* ref,
            Uint ref_len,
            double (*mismatch)(const void*, size_t, size_t),
            double (*delta)(bool));
    bool align_full_primal(bool force_align,
            Suffixtree* stree,
            Uchar* ref,
            Uint ref_len,
            double (*mismatch)(const void*, size_t, size_t),
            double (*delta)(bool));
    bool align_full_with_delta_primal(bool force_align,
            nucmer_Deltas& deltas);

    bool simple_align_left_primal(bool force_align,
            size_t rightmost_pos,
            Suffixtree* stree,
            Uchar* ref,
            Uint ref_len,
            double (*mismatch)(const void*, size_t, size_t),
            double (*delta)(bool));
    bool align_left_primal(bool force_align,
            size_t rightmost_pos,
            Suffixtree* stree,
            Uchar* ref,
            Uint ref_len,
            double (*mismatch)(const void*, size_t, size_t),
            double (*delta)(bool));
    bool align_left_with_delta_primal(bool force_align,
            size_t rightmost_pos,
            nucmer_Deltas& deltas);

    bool simple_align_right_primal(bool force_align,
            size_t leftmost_pos,
            Suffixtree* stree,
            Uchar* ref,
            Uint ref_len,
            double (*mismatch)(const void*, size_t, size_t),
            double (*delta)(bool));
    bool align_right_primal(bool force_align,
            size_t leftmost_pos,
            Suffixtree* stree,
            Uchar* ref,
            Uint ref_len,
            double (*mismatch)(const void*, size_t, size_t),
            double (*delta)(bool));
    bool align_right_with_delta_primal(bool force_align,
            size_t leftmost_pos,
            nucmer_Deltas& deltas);

    /*---------------- This is the dividing line --------------------*/

    long long compute_right_padding_len(const tigrinc::Nucmer_Cluster& clu);
    void merge_clusters(tigrinc::Nucmer_Cluster& clu1, tigrinc::Nucmer_Cluster& clu2,
            long long head_add, long long offset_1 = 0, long long offset_2 = 0);

    // res_indels must be empty when calling the following 3 func
    void nw_full(std::vector<long>& res_indels, 
            long long x_start, long long x_end, 
            long long y_start, long long y_end,
            long long head_add, long long& tail_add);
    void nw_left_unbound(std::vector<long>& res_indels,
            long long x_t, long long y_t,
            size_t& x_start, size_t& y_start,
            long long head_add, long long& tail_add); // head_add should always be 0 in this func call
    void nw_right_unbound(std::vector<long>& res_indels,
            long long x_s, long long y_s,
            size_t& x_end, size_t& y_end,
            long long head_add, long long& tail_add);

    static double default_mismatch(const void* obj, size_t px, size_t py);
    static double default_delta(bool is_insertion);
    static double max(const double aa, const double bb, const double cc);
};
}

#endif
