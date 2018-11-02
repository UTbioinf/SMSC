#include "Aligner.h"

namespace loon
{
//#define SAVE_NW

#ifdef SAVE_NW
static std::ofstream fout_save_nw("save_nw.log");
#endif

size_t Aligner::uabsminus(size_t a, size_t b)
{
    return ((a > b) ? (a - b) : (b - a));
}

Aligner::Aligner(size_t max_size/* =4096 */): 
        XX(NULL), YY(NULL), a(max_size),
        max_gap(20), max_diff(10),
        simple_align(false)
{
    cmp_obj = static_cast<const void*>(this);
}

void Aligner::set_simple_align()
{
    simple_align = true;
}

void Aligner::unset_simple_align()
{
    simple_align = false;
}

void Aligner::set_strs(const string& a, const string& b)
{
    XX = &a;
    YY = &b;
}

void Aligner::set_ref(const string& a)
{
    XX = &a;
}

void Aligner::set_qry(const string& b)
{
    YY = &b;
}

void Aligner::set_cmp_obj()
{
    cmp_obj = static_cast<const void*>(this);
}

void Aligner::set_cmp_obj(const void* obj)
{
    cmp_obj = obj;
}

void Aligner::set_parameters(size_t max_gap/* = 20*/, size_t max_diff/* = 10*/)
{
    this->max_gap = max_gap;
    this->max_diff = max_diff;
}

bool Aligner::fill_gaps(bool fill_ends, bool force_align, tigrinc::Nucmer_Delta& d,
        double (*mismatch)(const void*, size_t, size_t)/* = default_mismatch */, 
        double (*delta)(bool)/* = default_delta */)
{
    indels.clear();
    this->mismatch = mismatch;
    this->delta = delta;

    if(d.clusters.empty())
    {
        error_code = 1; // no alignment
        return false;
    }

    if(d.clusters.front().ref_start >= d.clusters.front().qry_start)
        error_code = 1<<2;
    else
        error_code = 0;
    if(XX->length() - d.clusters.back().ref_end >= YY->length() - d.clusters.back().qry_end)
        error_code |= 1<<1;    

    long long head_add = compute_right_padding_len( d.clusters[0] ), head_add_i;
    size_t cur_piece = 0;

    // fill gaps
    size_t largest_piece_id = 0;
    size_t largest_size = 0;
    long long largest_head_add = INF;
    for(size_t i = 1; i < d.clusters.size(); ++i)
    {
        head_add_i = compute_right_padding_len( d.clusters[i] );
        long long X_gap = d.clusters[i].ref_start - d.clusters[cur_piece].ref_end - 1;
        long long Y_gap = d.clusters[i].qry_start - d.clusters[cur_piece].qry_end - 1;
        if(abs(X_gap - Y_gap) <= max_diff && 
                (force_align || (!force_align && std::max(X_gap, Y_gap) <= max_gap) ))
        {
            merge_clusters(d.clusters[cur_piece], d.clusters[i], head_add);
            d.clusters[i].is_alive = false;
        }
        else
        {
            if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
            {
                largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
                largest_piece_id = cur_piece;
                largest_head_add = head_add;
            }
            cur_piece = i;
        }
        head_add = head_add_i;
    }
    if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
    {
        largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
        largest_piece_id = cur_piece;
        largest_head_add = head_add;
    }
    else if(largest_head_add == INF)
    {
        largest_head_add = head_add;
    }

    if(fill_ends)
    {
        // choose the largest alignment and fill the left and right padding
        long long left_padding = std::min(static_cast<long long>(d.clusters[ largest_piece_id ].ref_start), 
                static_cast<long long>(d.clusters[ largest_piece_id ].qry_start)) - 1;
        long long right_padding = std::min(static_cast<long long>(d.ref_len - d.clusters[ largest_piece_id ].ref_end), 
                static_cast<long long>(d.qry_len - d.clusters[ largest_piece_id ].qry_end));
        if(std::max(left_padding, right_padding) <= max_end_padding)
        {
            nw_left_unbound(indels, 
                    d.clusters[largest_piece_id].ref_start - 2, d.clusters[largest_piece_id].qry_start - 2,
                    XX_start, YY_start,
                    0, head_add);
            DEBUG_PRINT("XX_start = %lu, YY_start = %lu", XX_start, YY_start);
            ++XX_start, ++YY_start;

            if(!d.clusters[largest_piece_id].indels.empty())
            {
                if(d.clusters[largest_piece_id].indels[0] > 0)
                    d.clusters[largest_piece_id].indels[0] += head_add;
                else
                    d.clusters[largest_piece_id].indels[0] -= head_add;
            }
            indels.insert( indels.end(), d.clusters[largest_piece_id].indels.begin(), d.clusters[largest_piece_id].indels.end() );
            std::vector<long> res_indels;
            nw_right_unbound(res_indels,
                    d.clusters[largest_piece_id].ref_end, d.clusters[largest_piece_id].qry_end,
                    XX_end, YY_end,
                    largest_head_add, head_add);
            DEBUG_PRINT("XX_end = %lu, YY_end = %lu", XX_end, YY_end);
            ++XX_end, ++YY_end;

            indels.insert( indels.end(), res_indels.begin(), res_indels.end() );

            error_code = 0;
            return true;
        }
    }
    else
    {
        const tigrinc::Nucmer_Cluster& cc = d.clusters[ largest_piece_id ];
        indels.assign( cc.indels.begin(), cc.indels.end() );
        XX_start = cc.ref_start;    XX_end = cc.ref_end;
        YY_start = cc.qry_start;    YY_end = cc.qry_end;
        error_code = 0;
        return true;
    }
    return false;
}

bool Aligner::align_full(void (*compute_deltas_func)(const std::string&,
                const std::string&, nucmer_Deltas&, size_t, int),
        bool force_align/* = false */,
        double (*mismatch)(const void*, size_t, size_t)/* = default_mismatch */, 
        double (*delta)(bool)/* = default_delta */)
{
    indels.clear();
    this->mismatch = mismatch;
    this->delta = delta;

    nucmer_Deltas deltas;

    compute_deltas_func(*XX, *YY, deltas, 0, 3);
    return align_full_with_delta_primal(force_align, deltas);
}

bool Aligner::align_full(bool force_align/* = false*/,
        double (*mismatch)(const void*, size_t, size_t)/* = default_mismatch */, 
        double (*delta)(bool)/* = default_delta */)
{
    if(simple_align)
        return simple_align_full_primal(force_align, NULL, NULL, 0, mismatch, delta);
    else
        return align_full_primal(force_align, NULL, NULL, 0, mismatch, delta);
}

bool Aligner::align_full(Suffixtree& stree,
        Uchar* ref,
        Uint ref_len,
        bool force_align/* = false*/,
        double (*mismatch)(const void*, size_t, size_t)/* = default_mismatch */, 
        double (*delta)(bool)/* = default_delta */)
{
    if(simple_align)
        return simple_align_full_primal(force_align, &stree, ref, ref_len, mismatch, delta);
    else
        return align_full_primal(force_align, &stree, ref, ref_len, mismatch, delta);
}


bool Aligner::align_left(void (*compute_deltas_func)(const std::string&,
                const std::string&, nucmer_Deltas&, size_t, int),
        bool force_align/* = false */,
        size_t rightmost_pos/* = INF */,
        double (*mismatch)(const void*, size_t, size_t)/* = default_mismatch*/,
        double (*delta)(bool)/* = default_delta */)
{
    indels.clear();
    this->mismatch = mismatch;
    this->delta = delta;

    if(rightmost_pos == INF)
    {
        rightmost_pos = std::min(XX->length(), YY->length()) - 1;
    }

    nucmer_Deltas deltas;
    compute_deltas_func(*XX, *YY, deltas, rightmost_pos, 2);
    return align_left_with_delta_primal(force_align, rightmost_pos, deltas);
}

bool Aligner::align_left(bool force_align/* = false */,
        size_t rightmost_pos/* = INF */,
        double (*mismatch)(const void*, size_t, size_t)/* = default_mismatch*/,
        double (*delta)(bool)/* = default_delta */)
{
    if(simple_align)
        return simple_align_left_primal(force_align, rightmost_pos, NULL, NULL, 0, mismatch, delta);
    else
        return align_left_primal(force_align, rightmost_pos, NULL, NULL, 0, mismatch, delta);
}

bool Aligner::align_left(Suffixtree& stree,
        Uchar* ref,
        Uint ref_len,
        bool force_align/* = false*/,
        size_t rightmost_pos/* = INF*/,
        double (*mismatch)(const void*, size_t, size_t)/* = default_mismatch*/,
        double (*delta)(bool)/* = default_delta*/)
{
    ERROR_PRINT("Error: align_left() function with suffix tree provided is not implemented correctly!!!");
    exit(1);

    if(simple_align)
        return simple_align_left_primal(force_align, rightmost_pos, &stree, ref, ref_len, mismatch, delta);
    else
        return align_left_primal(force_align, rightmost_pos, &stree, ref, ref_len, mismatch, delta);
}

bool Aligner::align_right(void (*compute_deltas_func)(const std::string&,
                const std::string&, nucmer_Deltas&, size_t, int),
        bool force_align/* = false */,
        size_t leftmost_pos/* = INF */,
        double (*mismatch)(const void*, size_t, size_t)/* = default_mismatch*/,
        double (*delta)(bool)/* = defaul_delta*/)
{
    indels.clear();
    this->mismatch = mismatch;
    this->delta = delta;

    if(leftmost_pos == INF)
    {
        leftmost_pos = (XX->length() > YY->length() ? XX->length() - YY->length() : 0);
    }
    nucmer_Deltas deltas;
    compute_deltas_func(*XX, *YY, deltas, leftmost_pos, 1);

    return align_right_with_delta_primal(force_align, leftmost_pos, deltas);
}
        

bool Aligner::align_right(bool force_align/* = false*/,
        size_t leftmost_pos/* = INF*/,
        double (*mismatch)(const void*, size_t, size_t)/* = default_mismatch*/,
        double (*delta)(bool)/* = defaul_delta*/)
{
    if(simple_align)
        return simple_align_right_primal(force_align, leftmost_pos, NULL, NULL, 0, mismatch, delta);
    else
        return align_right_primal(force_align, leftmost_pos, NULL, NULL, 0, mismatch, delta);
}

bool Aligner::align_right(Suffixtree& stree,
        Uchar* ref,
        Uint ref_len,
        bool force_align/* = false */,
        size_t leftmost_pos/* = INF */,
        double (*mismatch)(const void*, size_t, size_t)/* = default_mismatch */,
        double (*delta)(bool)/* = default_delta */)
{
    ERROR_PRINT("Error: align_right() function with suffix tree provided is not implemented correctly!!!");
    exit(1);

    if(simple_align)
        return simple_align_right_primal(force_align, leftmost_pos, &stree, ref, ref_len, mismatch, delta);
    else
        return align_right_primal(force_align, leftmost_pos, &stree, ref, ref_len, mismatch, delta);
}

/*---------------- This is the dividing line --------------------*/

bool Aligner::simple_align_full_primal(bool force_align,
        Suffixtree* stree,
        Uchar* ref,
        Uint ref_len,
        double (*mismatch)(const void*, size_t, size_t), 
        double (*delta)(bool))
{
    this->mismatch = mismatch;
    this->delta = delta;
    const std::string& X = *XX;
    const std::string& Y = *YY;

    nucmer_Params nuc_p;
    nuc_p.minmatchlength = 5;
    nucmer_set_parameters( nuc_p );

    nucmer_Deltas deltas;
    if(stree == NULL)
        nucmer_alignall(X, Y, deltas);
    else
        nucmer_compute_deltas_once(*stree, ref, ref_len, Y, deltas, true);
    nucmer_reset_parameters();

    if(deltas.empty() || deltas[0].clusters.empty())
    {
        error_code = 1; // no alignment
        return false;
    }

    // merge 
    if(force_align && max_diff < default_force_max_diff)
        max_diff = default_force_max_diff;

    tigrinc::Nucmer_Delta& d = deltas[0];

    if(d.clusters.front().ref_start >= d.clusters.front().qry_start)
        error_code = 1<<2;
    else
        error_code = 0;
    if(XX->length() - d.clusters.back().ref_end >= YY->length() - d.clusters.back().qry_end)
        error_code |= 1<<1;

    size_t cur_piece = 0;

    // fill gaps
    size_t largest_piece_id = 0;
    size_t largest_size = 0;
    for(size_t i = 1; i < d.clusters.size(); ++i)
    {
        long long X_gap = d.clusters[i].ref_start - d.clusters[cur_piece].ref_end - 1;
        long long Y_gap = d.clusters[i].qry_start - d.clusters[cur_piece].qry_end - 1;
        if(abs(X_gap - Y_gap) <= max_diff && 
                (force_align || (!force_align && std::max(X_gap, Y_gap) <= max_gap) ))
        {
            d.clusters[cur_piece].ref_end = d.clusters[i].ref_end;
            d.clusters[cur_piece].qry_end = d.clusters[i].qry_end;
            d.clusters[i].is_alive = false;
        }
        else
        {
            if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
            {
                largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
                largest_piece_id = cur_piece;
            }
            cur_piece = i;
        }
    }
    if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
    {
        largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
        largest_piece_id = cur_piece;
    }

    // choose the largest alignment and fill the left and right padding
    long long left_padding = std::min(static_cast<long long>(d.clusters[ largest_piece_id ].ref_start), 
            static_cast<long long>(d.clusters[ largest_piece_id ].qry_start)) - 1;
    long long right_padding = std::min(static_cast<long long>(d.ref_len - d.clusters[ largest_piece_id ].ref_end), 
            static_cast<long long>(d.qry_len - d.clusters[ largest_piece_id ].qry_end));
    if(std::max(left_padding, right_padding) <= max_end_padding)
    {
        error_code = 0;
        return true;
    }

    return false;
}

bool Aligner::align_full_primal(bool force_align, 
        Suffixtree* stree,
        Uchar* ref,
        Uint ref_len,
        double (*mismatch)(const void*, size_t, size_t), 
        double (*delta)(bool))
{
    indels.clear();

    this->mismatch = mismatch;
    this->delta = delta;
    const std::string& X = *XX;
    const std::string& Y = *YY;

    nucmer_Params nuc_p;
    if(force_align)
    {
        nuc_p.minmatchlength = 5;
        nucmer_set_parameters( nuc_p );
    }
    else
    {
        nuc_p.minmatchlength = 20;
    }

    nucmer_Deltas deltas;
    if(stree == NULL)
        nucmer_alignall(X, Y, deltas);
    else
        nucmer_compute_deltas_once(*stree, ref, ref_len, Y, deltas, true);
    nucmer_reset_parameters();

    return align_full_with_delta_primal(force_align, deltas);
}

bool Aligner::align_full_with_delta_primal(bool force_align,
        nucmer_Deltas& deltas)
{
    if(deltas.empty() || deltas[0].clusters.empty())
    {
        error_code = 1; // no alignment
        return false;
    }

    // merge 
    if(force_align && max_diff < default_force_max_diff)
        max_diff = default_force_max_diff;

    tigrinc::Nucmer_Delta& d = deltas[0];

    if(d.clusters.front().ref_start >= d.clusters.front().qry_start)
        error_code = 1<<2;
    else
        error_code = 0;
    if(XX->length() - d.clusters.back().ref_end >= YY->length() - d.clusters.back().qry_end)
        error_code |= 1<<1;

    long long head_add = compute_right_padding_len( d.clusters[0] ), head_add_i;
    size_t cur_piece = 0;

    // fill gaps
    size_t largest_piece_id = 0;
    size_t largest_size = 0;
    long long largest_head_add = INF;
    for(size_t i = 1; i < d.clusters.size(); ++i)
    {
        head_add_i = compute_right_padding_len( d.clusters[i] );
        long long X_gap = d.clusters[i].ref_start - d.clusters[cur_piece].ref_end - 1;
        long long Y_gap = d.clusters[i].qry_start - d.clusters[cur_piece].qry_end - 1;
        if(abs(X_gap - Y_gap) <= max_diff && 
                (force_align || (!force_align && std::max(X_gap, Y_gap) <= max_gap) ))
        {
            merge_clusters(d.clusters[cur_piece], d.clusters[i], head_add);
            d.clusters[i].is_alive = false;
        }
        else
        {
            if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
            {
                largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
                largest_piece_id = cur_piece;
                largest_head_add = head_add;
            }
            cur_piece = i;
        }
        head_add = head_add_i;
    }
    if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
    {
        largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
        largest_piece_id = cur_piece;
        largest_head_add = head_add;
    }
    else if(largest_head_add == INF)
    {
        largest_head_add = head_add;
    }

    // choose the largest alignment and fill the left and right padding
    long long left_padding = std::min(static_cast<long long>(d.clusters[ largest_piece_id ].ref_start), 
            static_cast<long long>(d.clusters[ largest_piece_id ].qry_start)) - 1;
    long long right_padding = std::min(static_cast<long long>(d.ref_len - d.clusters[ largest_piece_id ].ref_end), 
            static_cast<long long>(d.qry_len - d.clusters[ largest_piece_id ].qry_end));
    if(std::max(left_padding, right_padding) <= max_end_padding)
    {
        nw_left_unbound(indels, 
                d.clusters[largest_piece_id].ref_start - 2, d.clusters[largest_piece_id].qry_start - 2,
                XX_start, YY_start,
                0, head_add);
        DEBUG_PRINT("XX_start = %lu, YY_start = %lu", XX_start, YY_start);
        ++XX_start, ++YY_start;

        if(!d.clusters[largest_piece_id].indels.empty())
        {
            if(d.clusters[largest_piece_id].indels[0] > 0)
                d.clusters[largest_piece_id].indels[0] += head_add;
            else
                d.clusters[largest_piece_id].indels[0] -= head_add;
        }
        indels.insert( indels.end(), d.clusters[largest_piece_id].indels.begin(), d.clusters[largest_piece_id].indels.end() );
        std::vector<long> res_indels;
        nw_right_unbound(res_indels,
                d.clusters[largest_piece_id].ref_end, d.clusters[largest_piece_id].qry_end,
                XX_end, YY_end,
                largest_head_add, head_add);
        DEBUG_PRINT("XX_end = %lu, YY_end = %lu", XX_end, YY_end);
        ++XX_end, ++YY_end;

        indels.insert( indels.end(), res_indels.begin(), res_indels.end() );

        error_code = 0;
        return true;
    }

    return false;
    
}

bool Aligner::simple_align_left_primal(bool force_align,
            size_t rightmost_pos,
            Suffixtree* stree,
            Uchar* ref,
            Uint ref_len,
            double (*mismatch)(const void*, size_t, size_t),
            double (*delta)(bool))
{
    indels.clear();

    this->mismatch = mismatch;
    this->delta = delta;

    if(rightmost_pos == INF)
    {
        rightmost_pos = std::min(XX->length(), YY->length()) - 1;
    }

    const std::string X = XX->substr(0, rightmost_pos + 1);
    const std::string& Y = *YY;

    nucmer_Params nuc_p;
    nuc_p.minmatchlength = 5;
    nucmer_set_parameters( nuc_p );

    nucmer_Deltas deltas;
    if(stree == NULL)
        nucmer_alignall(X, Y, deltas);
    else
        nucmer_compute_deltas_once(*stree, ref, ref_len, Y, deltas, true);
    nucmer_reset_parameters();

    if(deltas.empty() || deltas[0].clusters.empty())
    {
        error_code = 1;
        return false;
    }


    // merge 
    if(force_align && max_diff < default_force_max_diff)
        max_diff = default_force_max_diff;

    tigrinc::Nucmer_Delta& d = deltas[0];

    if(d.clusters.front().ref_start >= d.clusters.front().qry_start)
        error_code = 1<<2;
    else
        error_code = 0;
    if(XX->length() - d.clusters.back().ref_end >= YY->length() - d.clusters.back().qry_end)
        error_code |= 1<<1;

    size_t cur_piece = 0;

    // fill gaps
    size_t largest_piece_id = 0;
    size_t largest_size = 0;
    for(size_t i = 1; i < d.clusters.size(); ++i)
    {
        long long X_gap = d.clusters[i].ref_start - d.clusters[cur_piece].ref_end - 1;
        long long Y_gap = d.clusters[i].qry_start - d.clusters[cur_piece].qry_end - 1;
        if(abs(X_gap - Y_gap) <= max_diff && 
                (force_align || (!force_align && std::max(X_gap, Y_gap) <= max_gap) ))
        {
            d.clusters[cur_piece].ref_end = d.clusters[i].ref_end;
            d.clusters[cur_piece].qry_end = d.clusters[i].qry_end;
            d.clusters[i].is_alive = false;
        }
        else
        {
            if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
            {
                largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
                largest_piece_id = cur_piece;
            }
            cur_piece = i;
        }
    }
    if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
    {
        largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
        largest_piece_id = cur_piece;
    }

    // choose the largest alignment and fill the left and right padding
    long long left_padding = std::min(static_cast<long long>(d.clusters[ largest_piece_id ].ref_start), 
            static_cast<long long>(d.clusters[ largest_piece_id ].qry_start)) - 1;
    long long right_padding = std::min(static_cast<long long>(XX->length() - d.clusters[ largest_piece_id ].ref_end), 
            static_cast<long long>(d.qry_len - d.clusters[ largest_piece_id ].qry_end));
    if(std::max(left_padding, right_padding) <= max_end_padding)
    {
        return true;
    }
    return false;
}

bool Aligner::align_left_primal(bool force_align,
            size_t rightmost_pos,
            Suffixtree* stree,
            Uchar* ref,
            Uint ref_len,
            double (*mismatch)(const void*, size_t, size_t),
            double (*delta)(bool))
{
    indels.clear();

    this->mismatch = mismatch;
    this->delta = delta;

    if(rightmost_pos == INF)
    {
        rightmost_pos = std::min(XX->length(), YY->length()) - 1;
    }

    const std::string X = XX->substr(0, rightmost_pos + 1);
    const std::string& Y = *YY;

    nucmer_Params nuc_p;
    nuc_p.minmatchlength = 5;
    nucmer_set_parameters( nuc_p );

    nucmer_Deltas deltas;
    if(stree == NULL)
        nucmer_alignall(X, Y, deltas);
    else
        nucmer_compute_deltas_once(*stree, ref, ref_len, Y, deltas, true);
    nucmer_reset_parameters();

    return align_left_with_delta_primal(force_align, rightmost_pos, deltas);

}

bool Aligner::align_left_with_delta_primal(bool force_align,
        size_t rightmost_pos,
        nucmer_Deltas& deltas)
{
    if(deltas.empty() || deltas[0].clusters.empty())
    {
        error_code = 1;
        return false;
    }


    // merge 
    if(force_align && max_diff < default_force_max_diff)
        max_diff = default_force_max_diff;

    tigrinc::Nucmer_Delta& d = deltas[0];

    if(d.clusters.front().ref_start >= d.clusters.front().qry_start)
        error_code = 1<<2;
    else
        error_code = 0;
    if(XX->length() - d.clusters.back().ref_end >= YY->length() - d.clusters.back().qry_end)
        error_code |= 1<<1;

    long long head_add = compute_right_padding_len( d.clusters[0] ), head_add_i;
    size_t cur_piece = 0;

    // fill gaps
    size_t largest_piece_id = 0;
    size_t largest_size = 0;
    long long largest_head_add = INF;
    for(size_t i = 1; i < d.clusters.size(); ++i)
    {
        head_add_i = compute_right_padding_len( d.clusters[i] );
        long long X_gap = d.clusters[i].ref_start - d.clusters[cur_piece].ref_end - 1;
        long long Y_gap = d.clusters[i].qry_start - d.clusters[cur_piece].qry_end - 1;
        if(abs(X_gap - Y_gap) <= max_diff && 
                (force_align || (!force_align && std::max(X_gap, Y_gap) <= max_gap) ))
        {
            merge_clusters(d.clusters[cur_piece], d.clusters[i], head_add);
            d.clusters[i].is_alive = false;
        }
        else
        {
            if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
            {
                largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
                largest_piece_id = cur_piece;
                largest_head_add = head_add;
            }
            cur_piece = i;
        }
        head_add = head_add_i;
    }
    if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
    {
        largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
        largest_piece_id = cur_piece;
        largest_head_add = head_add;
    }
    else if(largest_head_add == INF)
    {
        largest_head_add = head_add;
    }

    // choose the largest alignment and fill the left and right padding
    long long left_padding = std::min(static_cast<long long>(d.clusters[ largest_piece_id ].ref_start), 
            static_cast<long long>(d.clusters[ largest_piece_id ].qry_start)) - 1;
    long long right_padding = std::min(static_cast<long long>(XX->length() - d.clusters[ largest_piece_id ].ref_end), 
            static_cast<long long>(d.qry_len - d.clusters[ largest_piece_id ].qry_end));
    if(std::max(left_padding, right_padding) <= max_end_padding)
    {
        nw_left_unbound(indels, 
                d.clusters[largest_piece_id].ref_start - 2, 
                d.clusters[largest_piece_id].qry_start - 2,
                XX_start, YY_start,
                0, head_add);
        DEBUG_PRINT("XX_start = %lu, YY_start = %lu", XX_start, YY_start);
        ++XX_start, ++YY_start;

        if(!d.clusters[largest_piece_id].indels.empty())
        {
            if(d.clusters[largest_piece_id].indels[0] > 0)
                d.clusters[largest_piece_id].indels[0] += head_add;
            else
                d.clusters[largest_piece_id].indels[0] -= head_add;
        }
        indels.insert( indels.end(), d.clusters[largest_piece_id].indels.begin(), d.clusters[largest_piece_id].indels.end() );
        std::vector<long> res_indels;
        nw_right_unbound(res_indels,
                d.clusters[largest_piece_id].ref_end, 
                d.clusters[largest_piece_id].qry_end,
                XX_end, YY_end,
                largest_head_add, head_add);
        DEBUG_PRINT("XX_end = %lu, YY_end = %lu", XX_end, YY_end);
        ++XX_end, ++YY_end;

        indels.insert( indels.end(), res_indels.begin(), res_indels.end() );

        return true;
    }
    return false;
}

bool Aligner::simple_align_right_primal(bool force_align,
        size_t leftmost_pos,
        Suffixtree* stree,
        Uchar* ref,
        Uint ref_len,
        double (*mismatch)(const void*, size_t, size_t),
        double (*delta)(bool))
{
    indels.clear();
    
    this->mismatch = mismatch;
    this->delta = delta;

    if(leftmost_pos == INF)
    {
        leftmost_pos = (XX->length() > YY->length() ? XX->length() - YY->length() : 0);
    }

    const std::string X = XX->substr(leftmost_pos);
    //const std::string& Y = *YY;
    const std::string& Y = YY->length() > X.length() ? YY->substr(0, X.length()) : *YY;

    if(X.length() <= 0 || Y.length() <=0)
    {
        error_code = ((!X.length() << 1) | (!Y.length())) << 3;
        return false;
    }

    nucmer_Params nuc_p;
    nuc_p.minmatchlength = 5;
    nucmer_set_parameters( nuc_p );

    nucmer_Deltas deltas;
    if(stree == NULL)
        nucmer_alignall(X, Y, deltas);
    else
        nucmer_compute_deltas_once(*stree, ref, ref_len, Y, deltas, true);
    nucmer_reset_parameters();

    if(deltas.empty() || deltas[0].clusters.empty())
    {
        error_code = 1;
        return false;
    }

    // merge 
    if(force_align && max_diff < default_force_max_diff)
        max_diff = default_force_max_diff;

    tigrinc::Nucmer_Delta& d = deltas[0];

    if(d.clusters.front().ref_start + leftmost_pos >= d.clusters.front().qry_start)
        error_code = 1<<2;
    else
        error_code = 0;
    if(XX->length() - d.clusters.back().ref_end - leftmost_pos >= YY->length() - d.clusters.back().qry_end)
        error_code |= 1<<1;

    size_t cur_piece = 0;

    // fill gaps
    const size_t piece_size_threshold = 10;
    const long long INFLongLong = std::numeric_limits<long long>::max(); 
    long long min_padding = INFLongLong;
    size_t largest_piece_id = 0;
    size_t largest_size = 0;
    long long left_padding, right_padding;
    for(size_t i = 1; i < d.clusters.size(); ++i)
    {
        long long X_gap = d.clusters[i].ref_start - d.clusters[cur_piece].ref_end - 1;
        long long Y_gap = d.clusters[i].qry_start - d.clusters[cur_piece].qry_end - 1;
        if(abs(X_gap - Y_gap) <= max_diff && 
                (force_align || (!force_align && std::max(X_gap, Y_gap) <= max_gap) ))
        {
            d.clusters[cur_piece].ref_end = d.clusters[i].ref_end;
            d.clusters[cur_piece].qry_end = d.clusters[i].qry_end;
            d.clusters[i].is_alive = false;
        }
        else
        {
            left_padding = std::min(static_cast<long long>(d.clusters[ cur_piece ].ref_start + leftmost_pos), 
                    static_cast<long long>(d.clusters[ cur_piece ].qry_start) - 1);
            right_padding = std::min(static_cast<long long>(d.ref_len - d.clusters[ cur_piece ].ref_end), 
                    static_cast<long long>(d.qry_len - d.clusters[ cur_piece ].qry_end));
            //if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
            if(min_padding > std::max(left_padding, right_padding))
            {
                min_padding = std::max(left_padding, right_padding);
                largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
                largest_piece_id = cur_piece;
            }
            cur_piece = i;
        }
    }
    left_padding = std::min(static_cast<long long>(d.clusters[ cur_piece ].ref_start + leftmost_pos), 
            static_cast<long long>(d.clusters[ cur_piece ].qry_start) - 1);
    right_padding = std::min(static_cast<long long>(d.ref_len - d.clusters[ cur_piece ].ref_end), 
            static_cast<long long>(d.qry_len - d.clusters[ cur_piece ].qry_end));
    //if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
    if(min_padding > std::max(left_padding, right_padding))
    {
        min_padding = std::max(left_padding, right_padding);
        largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
        largest_piece_id = cur_piece;
    }

    // choose the largest alignment and fill the left and right padding
    // long long left_padding = std::min(static_cast<long long>(d.clusters[ largest_piece_id ].ref_start + leftmost_pos), 
    //         static_cast<long long>(d.clusters[ largest_piece_id ].qry_start) - 1);
    // long long right_padding = std::min(static_cast<long long>(d.ref_len - d.clusters[ largest_piece_id ].ref_end), 
    //         static_cast<long long>(d.qry_len - d.clusters[ largest_piece_id ].qry_end));
    // if(std::max(left_padding, right_padding) <= max_end_padding)
    if(min_padding <= max_end_padding)
    {
        return true;
    }
    return false;
}

bool Aligner::align_right_primal(bool force_align,
        size_t leftmost_pos,
        Suffixtree* stree,
        Uchar* ref,
        Uint ref_len,
        double (*mismatch)(const void*, size_t, size_t),
        double (*delta)(bool))
{
    indels.clear();
    this->mismatch = mismatch;
    this->delta = delta;

    if(leftmost_pos == INF)
    {
        leftmost_pos = (XX->length() > YY->length() ? XX->length() - YY->length() : 0);
    }

    const std::string X = XX->substr(leftmost_pos);
    const std::string& Y = YY->length() > X.length() ? YY->substr(0, X.length()) : *YY;

    if(X.length() <= 0 || Y.length() <=0)
    {
        error_code = ((!X.length() << 1) | (!Y.length())) << 3;
        return false;
    }
    nucmer_Params nuc_p;
    nuc_p.minmatchlength = 5;
    nucmer_set_parameters( nuc_p );

    nucmer_Deltas deltas;
    if(stree == NULL)
        nucmer_alignall(X, Y, deltas);
    else
        nucmer_compute_deltas_once(*stree, ref, ref_len, Y, deltas, true);
    nucmer_reset_parameters();

    
    return align_right_with_delta_primal(force_align, leftmost_pos, deltas);

}

bool Aligner::align_right_with_delta_primal(bool force_align,
        size_t leftmost_pos,
        nucmer_Deltas& deltas)
{
    if(deltas.empty() || deltas[0].clusters.empty())
    {
        error_code = 1;
        return false;
    }

#ifdef DEBUG
    deltas.check_deltas();
#endif
    // merge 
    if(force_align && max_diff < default_force_max_diff)
        max_diff = default_force_max_diff;

    tigrinc::Nucmer_Delta& d = deltas[0];

    if(d.clusters.front().ref_start + leftmost_pos >= d.clusters.front().qry_start)
        error_code = 1<<2;
    else
        error_code = 0;
    if(XX->length() - d.clusters.back().ref_end - leftmost_pos >= YY->length() - d.clusters.back().qry_end)
        error_code |= 1<<1;

    long long head_add = compute_right_padding_len( d.clusters[0] ), head_add_i;
    size_t cur_piece = 0;

    // fill gaps
    const size_t piece_size_threshold = 10;
    const long long INFLongLong = std::numeric_limits<long long>::max(); 
    long long min_padding = INFLongLong;
    size_t largest_piece_id = 0;
    size_t largest_size = 0;
    long long largest_head_add = INFLongLong;
    long long left_padding, right_padding;
    for(size_t i = 1; i < d.clusters.size(); ++i)
    {
        head_add_i = compute_right_padding_len( d.clusters[i] );
        long long X_gap = d.clusters[i].ref_start - d.clusters[cur_piece].ref_end - 1;
        long long Y_gap = d.clusters[i].qry_start - d.clusters[cur_piece].qry_end - 1;
        if(abs(X_gap - Y_gap) <= max_diff && 
                (force_align || (!force_align && std::max(X_gap, Y_gap) <= max_gap) ))
        {
            //DEBUG_PRINT("merge them");
            merge_clusters(d.clusters[cur_piece], d.clusters[i], head_add, leftmost_pos);
            d.clusters[i].is_alive = false;
        }
        else
        {
            left_padding = std::min(static_cast<long long>(d.clusters[ cur_piece ].ref_start + leftmost_pos), 
                    static_cast<long long>(d.clusters[ cur_piece ].qry_start) - 1);
            right_padding = std::min(static_cast<long long>(d.ref_len - d.clusters[ cur_piece ].ref_end), 
                    static_cast<long long>(d.qry_len - d.clusters[ cur_piece ].qry_end));
            //DEBUG_PRINT("max_diff = %lld", abs(X_gap - Y_gap));
            //if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
            //if((uabsminus(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1, largest_size) < piece_size_threshold && min_padding > std::max(left_padding, right_padding)) || (d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size))
            if(min_padding > std::max(left_padding, right_padding))
            {
                min_padding = std::max(left_padding, right_padding);
                largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
                largest_piece_id = cur_piece;
                largest_head_add = head_add;
            }
            cur_piece = i;
        }
        head_add = head_add_i;
    }
    left_padding = std::min(static_cast<long long>(d.clusters[ cur_piece ].ref_start + leftmost_pos), 
            static_cast<long long>(d.clusters[ cur_piece ].qry_start) - 1);
    right_padding = std::min(static_cast<long long>(d.ref_len - d.clusters[ cur_piece ].ref_end), 
            static_cast<long long>(d.qry_len - d.clusters[ cur_piece ].qry_end));
    //DEBUG_PRINT("max_diff = %lld", abs(X_gap - Y_gap));
    //if(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size)
    //if((uabsminus(d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1, largest_size) < piece_size_threshold && min_padding > std::max(left_padding, right_padding)) || (d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1 > largest_size))
    if(min_padding > std::max(left_padding, right_padding))
    {
        min_padding = std::max(left_padding, right_padding);
        largest_size = d.clusters[cur_piece].qry_end - d.clusters[cur_piece].qry_start + 1;
        largest_piece_id = cur_piece;
        largest_head_add = head_add;
    }
    else if(largest_head_add == INFLongLong)
    {
        largest_head_add = head_add;
    }

    // choose the largest alignment and fill the left and right padding
    //left_padding = std::min(static_cast<long long>(d.clusters[ largest_piece_id ].ref_start + leftmost_pos), 
    //         static_cast<long long>(d.clusters[ largest_piece_id ].qry_start) - 1);
    //right_padding = std::min(static_cast<long long>(d.ref_len - d.clusters[ largest_piece_id ].ref_end), 
    //         static_cast<long long>(d.qry_len - d.clusters[ largest_piece_id ].qry_end));
    //min_padding = std::max(left_padding, right_padding);
    //if(std::max(left_padding, right_padding) <= max_end_padding)
    if(min_padding <= max_end_padding)
    {
        nw_left_unbound(indels, 
                d.clusters[largest_piece_id].ref_start - 2 + leftmost_pos,
                d.clusters[largest_piece_id].qry_start - 2,
                XX_start, YY_start,
                0, head_add);
        DEBUG_PRINT("XX_start = %lu, YY_start = %lu", XX_start, YY_start);
        ++XX_start, ++YY_start;

        if(!d.clusters[largest_piece_id].indels.empty())
        {
            if(d.clusters[largest_piece_id].indels[0] > 0)
                d.clusters[largest_piece_id].indels[0] += head_add;
            else
                d.clusters[largest_piece_id].indels[0] -= head_add;
        }
        indels.insert( indels.end(), d.clusters[largest_piece_id].indels.begin(), d.clusters[largest_piece_id].indels.end() );
        std::vector<long> res_indels;
        nw_right_unbound(res_indels,
                d.clusters[largest_piece_id].ref_end + leftmost_pos,
                d.clusters[largest_piece_id].qry_end,
                XX_end, YY_end,
                largest_head_add, head_add);
        DEBUG_PRINT("XX_end = %lu, YY_end = %lu", XX_end, YY_end);
        ++XX_end, ++YY_end;

        indels.insert( indels.end(), res_indels.begin(), res_indels.end() );

        return true;
    }
    //ERROR_PRINT("left padding = %lld, right padding = %lld, min_padding = %lld, left_most = %lu, ref_left = %lld, qry_left = %lld", left_padding, right_padding, min_padding, leftmost_pos, d.clusters[ largest_piece_id ].ref_start + leftmost_pos, d.clusters[ largest_piece_id ].qry_start - 1);
    return false;
}

/*---------------- This is the dividing line --------------------*/

long long Aligner::compute_right_padding_len(const tigrinc::Nucmer_Cluster& clu)
{
    long ref_pos = -1;
    long qry_pos = -1;

    for(std::vector<long>::const_iterator indel_it = clu.indels.begin(); indel_it != clu.indels.end(); ++indel_it)
    {
        if(*indel_it > 0)
        {
            ref_pos += *indel_it;
            qry_pos += *indel_it - 1;
        }
        else
        {
            ref_pos -= *indel_it + 1;
            qry_pos -= *indel_it;
        }
    }
    //DEBUG_PRINT("clu_ref: (start, end, pos) = (%lld, %lld, %ld)", clu.ref_start, clu.ref_end, ref_pos);
    //DEBUG_PRINT("clu_qry: (start, end, pos) = (%lld, %lld, %ld)", clu.qry_start, clu.qry_end, qry_pos);
    assert(clu.ref_end - clu.ref_start - ref_pos == clu.qry_end - clu.qry_start - qry_pos);
    assert(clu.ref_end - clu.ref_start - ref_pos >= 0);
    return clu.ref_end - clu.ref_start - ref_pos;
}

void Aligner::merge_clusters(tigrinc::Nucmer_Cluster& clu1, tigrinc::Nucmer_Cluster& clu2, 
        long long head_add, long long offset_1/* = 0 */, long long offset_2/* = 0 */)
{
    std::vector<long> res_indels;
    nw_full(res_indels,
            clu1.ref_end + offset_1, clu2.ref_start - 2 + offset_1,
            clu1.qry_end + offset_2, clu2.qry_start - 2 + offset_2,
            head_add, head_add);

    clu1.indels.insert( clu1.indels.end(), res_indels.begin(), res_indels.end() );
    if(!clu2.indels.empty())
    {
        if(clu2.indels[0] > 0)
            clu2.indels[0] += head_add;
        else
            clu2.indels[0] -= head_add;
        clu1.indels.insert( clu1.indels.end(), clu2.indels.begin(), clu2.indels.end() );
    }

    clu1.ref_end = clu2.ref_end;
    clu1.qry_end = clu2.qry_end;
}

void Aligner::nw_full(std::vector<long>& res_indels, 
        long long x_start, long long x_end, 
        long long y_start, long long y_end,
        long long head_add, long long& tail_add)
{
    DEBUG_PRINT("call nw_full()");
    const std::string& X = *XX;
    const std::string& Y = *YY;
    long long m = x_end - x_start + 1;
    long long n = y_end - y_start + 1;

#ifdef SAVE_NW
    if(m > 0 || n > 0)
    {
        fout_save_nw << ">full " << m << " " << n << std::endl;
        fout_save_nw << X.substr(x_start, m) << std::endl;
        fout_save_nw << Y.substr(y_start, n) << std::endl;
    }
#endif
    
    if(m == 0)
    {
        if(n == 0)
            return;
        // deletion
        res_indels.resize(n, -1);
        res_indels[0] -= head_add;
        tail_add = 0;
        return;
    }
    else if(n == 0)
    {
        // insertion
        res_indels.resize(m, 1);
        res_indels[0] += head_add;
        tail_add = 0;
        return;
    }

    // run the algorithm
    a.resize( m+1, n+1 );
    for(long long i = 0; i <= m; ++i)
        a.at(i, 0) = i * delta(true); // insertion
    for(long long j = 0; j <= n; ++j)
        a.at(0, j) = j * delta(false); // deletion

    for(long long j = 1; j <= n; ++j)
        for(long long i = 1; i <= m; ++i)
        {
            a.at(i, j) = max(
                        a.at(i-1, j-1) + mismatch(cmp_obj, x_start + i - 1, y_start + j - 1),
                        a.at(i-1, j) + delta(true), // insertion
                        a.at(i, j-1) + delta(false) // deletion
                    );
        }

    // compute the deltas
    bool tail_add_done = false;
    bool last_is_insert;
    long indel_cnt = 1;
    while(m > 0 || n > 0)
    {
        if(m > 0 && n > 0 && 
                fabs(a.at(m, n) - a.at(m-1, n-1) - mismatch(cmp_obj, x_start + m - 1, y_start + n - 1)) < ZERO)
        {
            --m;
            --n;
            ++indel_cnt;
        }
        else if(m > 0 && fabs(a.at(m, n) - a.at(m-1, n) - delta(true)) <= ZERO)
        {
            --m;
            if(!tail_add_done)
            {
                tail_add_done = true;
                tail_add = indel_cnt - 1;
            }
            else
            {
                res_indels.push_back( last_is_insert ? indel_cnt : -indel_cnt);
            }

            last_is_insert = true;
            indel_cnt = 1;
        }
        else
        {
            --n;
            if(!tail_add_done)
            {
                tail_add_done = true;
                tail_add = indel_cnt - 1;
            }
            else
            {
                res_indels.push_back( last_is_insert ? indel_cnt : -indel_cnt );
            }
            last_is_insert = false;
            indel_cnt = 1;
        }
    }

    if(!tail_add_done)
    {
        tail_add = indel_cnt - 1 + head_add;
    }
    else
    {
        res_indels.push_back( last_is_insert ? indel_cnt + head_add : -indel_cnt - head_add );
        reverse(res_indels.begin(), res_indels.end());
    }
}

void Aligner::nw_left_unbound(std::vector<long>& res_indels,
        long long x_t, long long y_t,
        size_t& x_start, size_t& y_start,
        long long head_add, long long& tail_add)
{
    long long x_len = x_t + 1;
    long long y_len = y_t + 1;
    long long x_n, y_n;

    DEBUG_PRINT("nw_left_unbound: x_t = %lld, y_t = %lld", x_t, y_t);

#ifdef SAVE_NW
    if(x_len > 0 && y_len > 0)
    {
        fout_save_nw << ">left " << x_len << " " << y_len << std::endl;
        fout_save_nw << XX->substr(0, x_len) << std::endl;
        fout_save_nw << YY->substr(0, y_len) << std::endl;
    }
#endif
    if(x_len <= 0)
    {// no insertion or deletion required in this case
        x_start = 0;
        y_start = std::max(static_cast<size_t>(0), static_cast<size_t>(y_len));
        tail_add = 0;
        return;
    }
    else if(y_len <= 0)
    {
        x_start = std::max(static_cast<size_t>(0), static_cast<size_t>(x_len));
        y_start = 0;
        tail_add = 0;
        return;
    }

    size_t best_row = 0;
    if(x_len <= y_len)
    {
        a.resize_column(x_len + 1);

        DEBUG_PRINT("x_len <= y_len: init");
        for(size_t j = 0; j <= x_len; ++j)
            a.at(0, j) = j * delta(true);

        DEBUG_PRINT("x_len <= y_len: run algorithm");
        size_t smaller_cnt = 0;
        for(size_t i = 1; i <= y_len; ++i)
        {
            a.add_row();
            a.at(i, 0) = i * delta(false);

            for(size_t j=1; j <= x_len; ++j)
            {
                a.at(i, j) = max(
                            a.at(i-1, j-1) + mismatch(cmp_obj, x_len - j, y_len - i),
                            a.at(i-1, j) + delta(false),
                            a.at(i, j-1) + delta(true)
                        );
            }

            if(a.at(best_row, x_len) < a.at(i, x_len))
            {
                best_row = i;
                smaller_cnt = 0;
            }
            else if(++smaller_cnt > max_ext && i>x_len)
            {
                break;
            }
        }

        x_n = x_len;
        y_n = best_row;

        x_start = 0;
        y_start = y_len - best_row;

        DEBUG_PRINT("x_len <= y_len: compute indels");
        long indel_cnt = 1;
        while(x_n > 0 || y_n > 0)
        {
            if(y_n > 0 && x_n > 0 && fabs(a.at(y_n, x_n) - a.at(y_n - 1, x_n - 1) - mismatch(cmp_obj, x_len - x_n, y_len - y_n)) <= ZERO)
            {
                --y_n;
                --x_n;
                ++indel_cnt;
            }
            else if(y_n > 0 && fabs(a.at(y_n, x_n) - a.at(y_n - 1, x_n) - delta(false) ) <= ZERO)
            {
                --y_n;
                res_indels.push_back( -indel_cnt );
                indel_cnt = 1;
            }
            else
            {
                --x_n;
                res_indels.push_back( indel_cnt );
                indel_cnt = 1;
            }
        }
        if(res_indels.empty())
        {
            tail_add = indel_cnt - 1 + head_add;
        }
        else
        {
            res_indels[0] += head_add;
            tail_add = indel_cnt - 1;
        }
    }
    else
    {
        a.resize_column( y_len + 1 );

        for(size_t j=0; j <= y_len; ++j)
            a.at(0, j) = j * delta(false);

        size_t smaller_cnt = 0;
        for(size_t i = 1; i <= x_len; ++i)
        {
            a.add_row();
            a.at(i, 0) = i * delta(true);
            for(size_t j = 1; j <= y_len; ++j)
            {
                a.at(i, j) = max(
                            a.at(i-1, j-1) + mismatch(cmp_obj, x_len - i, y_len - j),
                            a.at(i-1, j) + delta(true),
                            a.at(i, j-1) + delta(false)
                        );
            }

            if(a.at(best_row, y_len) < a.at(i, y_len))
            {
                best_row = i;
                smaller_cnt = 0;
            }
            else if(++smaller_cnt > max_ext && i > y_len)
            {
                break;
            }
        }

        x_n = best_row;
        y_n = y_len;

        x_start = x_len - best_row;
        y_start = 0;


        long indel_cnt = 1;
        while(x_n > 0 || y_n > 0)
        {
            if(x_n > 0 && y_n > 0 && fabs(a.at(x_n, y_n) - a.at(x_n - 1, y_n - 1) - mismatch(cmp_obj, x_len - x_n, y_len - y_n)) <= ZERO)
            {
                --x_n;
                --y_n;
                ++indel_cnt;
            }
            else if(x_n > 0 && fabs(a.at(x_n, y_n) - a.at(x_n - 1, y_n) - delta(true)) <= ZERO)
            {
                --x_n;
                res_indels.push_back( indel_cnt );
                indel_cnt = 1;
            }
            else
            {
                --y_n;
                res_indels.push_back( -indel_cnt );
                indel_cnt = 1;
            }
        }

        if(res_indels.empty())
        {
            tail_add = indel_cnt - 1 + head_add;
        }
        else
        {
            res_indels[0] += head_add;
            tail_add = indel_cnt - 1;
        }
    }
}

void Aligner::nw_right_unbound(std::vector<long>& res_indels,
        long long x_s, long long y_s,
        size_t& x_end, size_t& y_end,
        long long head_add, long long& tail_add)
{
    long long x_len = XX->length() - x_s;
    long long y_len = YY->length() - y_s;

    long long x_n, y_n;

    DEBUG_PRINT("nw_right_unbound: x_s = %lld, y_s = %lld", x_s, y_s);

#ifdef SAVE_NW
    if(x_len > 0 && y_len > 0)
    {
        fout_save_nw << ">right " << x_len << " " << y_len << std::endl;
        fout_save_nw << XX->substr(x_s, x_len) << std::endl;
        fout_save_nw << YY->substr(y_s, y_len) << std::endl;
    }
#endif

    if(x_len <= 0 || y_len <= 0)
    {// no insertion or deletion required in this case
        x_end = x_s - 1;
        y_end = y_s - 1;
        tail_add = head_add;
        return;
    }

    size_t best_row = 0;

    if(x_len <= y_len)
    {
        a.resize_column( x_len + 1 );

        for(size_t j = 0; j <= x_len; ++j)
            a.at(0, j) = j * delta(true);

        size_t smaller_cnt = 0;
        for(size_t i = 1; i <= y_len; ++i)
        {
            a.add_row();
            a.at(i, 0) = i * delta(false);

            for(size_t j = 1; j <= x_len; ++j)
            {
                a.at(i, j) = max(
                            a.at(i-1, j-1) + mismatch(cmp_obj, x_s + j - 1, y_s + i - 1),
                            a.at(i-1, j) + delta(false),
                            a.at(i, j-1) + delta(true)
                        );
            }

            if(a.at(best_row, x_len) < a.at(i, x_len))
            {
                best_row = i;
                smaller_cnt = 0;
            }
            else if(++smaller_cnt > max_ext && i > x_len)
            {
                break;
            }
        }

        // compute deltas & (x_end, y_end)
        x_n = x_len;
        y_n = best_row;

        x_end = XX->length()- 1;
        y_end = y_s + best_row - 1;

        bool tail_add_done = false;
        bool last_is_insert;
        long indel_cnt = 1;
        while(y_n > 0 || x_n > 0)
        {
            if(y_n > 0 && x_n > 0 && fabs(a.at(y_n, x_n) - a.at(y_n - 1, x_n -1 ) - mismatch(cmp_obj, x_s + x_n - 1, y_s + y_n - 1)) <= ZERO)
            {
                --y_n;
                --x_n;
                ++indel_cnt;
            }
            else if(y_n > 0 && fabs(a.at(y_n, x_n) - a.at(y_n - 1, x_n) - delta(false)) <= ZERO)
            {
                --y_n;
                if(!tail_add_done)
                {
                    tail_add_done = true;
                    tail_add = indel_cnt - 1;
                }
                else
                {
                    res_indels.push_back( last_is_insert ? indel_cnt : -indel_cnt );
                }
                last_is_insert = false;
                indel_cnt = 1;
            }
            else
            {
                --x_n;
                if(!tail_add_done)
                {
                    tail_add_done = true;
                    tail_add = indel_cnt - 1;
                }
                else
                {
                    res_indels.push_back( last_is_insert ? indel_cnt : -indel_cnt );
                }
                last_is_insert = true;
                indel_cnt = 1;
            }
        }

        if(!tail_add_done)
        {
            tail_add = indel_cnt - 1 + head_add;
        }
        else
        {
            res_indels.push_back( last_is_insert ? indel_cnt + head_add : -indel_cnt - head_add );
            reverse(res_indels.begin(), res_indels.end());
        }
    }
    else
    {
        a.resize_column( y_len + 1 );

        for(size_t j = 0; j <= y_len; ++j)
            a.at(0, j) = j * delta(false);

        size_t smaller_cnt = 0;
        for(size_t i = 1; i <= x_len; ++i)
        {
            a.add_row();
            a.at(i, 0) = i * delta(true);
            for(size_t j = 1; j <= y_len; ++j)
            {
                a.at(i, j) = max(
                            a.at(i-1, j-1) + mismatch(cmp_obj, x_s + i - 1, y_s + j - 1),
                            a.at(i-1, j) + delta(true),
                            a.at(i, j-1) + delta(false)
                        );
            }

            if(a.at(best_row, y_len) < a.at(i, y_len))
            {
                best_row = i;
                smaller_cnt = 0;
            }
            else if(++smaller_cnt > max_ext && i > y_len)
            {
                break;
            }
        }

        // compute deltas & (x_end, y_end)
        x_n = best_row;
        y_n = y_len;

        x_end = x_s + best_row - 1;
        y_end = YY->length() - 1;

        bool tail_add_done = false;
        bool last_is_insert;
        long indel_cnt = 1;
        while(x_n > 0 || y_n > 0)
        {
            if(x_n > 0 && y_n > 0 && fabs(a.at(x_n, y_n) - a.at(x_n - 1, y_n - 1) - mismatch(cmp_obj, x_s + x_n - 1, y_s + y_n - 1)) <= ZERO)
            {
                --x_n;
                --y_n;
                ++indel_cnt;
            }
            else if(x_n > 0 && fabs(a.at(x_n, y_n) - a.at(x_n - 1, y_n) - delta(true)) <= ZERO)
            {
                --x_n;
                if(!tail_add_done)
                {
                    tail_add_done = true;
                    tail_add = indel_cnt - 1;
                }
                else
                {
                    res_indels.push_back( last_is_insert ? indel_cnt : -indel_cnt );
                }
                last_is_insert = true;
                indel_cnt = 1;
                
            }
            else
            {
                --y_n;
                if(!tail_add_done)
                {
                    tail_add_done = true;
                    tail_add = indel_cnt - 1;
                }
                else
                {
                    res_indels.push_back( last_is_insert ? indel_cnt : -indel_cnt );
                }

                last_is_insert = false;
                indel_cnt = 1;
            }
        }

        if(!tail_add_done)
        {
            tail_add = indel_cnt - 1 + head_add;
        }
        else
        {
            res_indels.push_back( last_is_insert ? indel_cnt + head_add : -indel_cnt - head_add );
            reverse(res_indels.begin(), res_indels.end());
        }
    }
}

/*static*/ double Aligner::default_mismatch(const void* obj, size_t px, size_t py)
{
    const Aligner* t_obj = static_cast<const Aligner*>(obj);
    if(t_obj->XX->at(px) == t_obj->YY->at(py))      return 3;
    else    return -2;
}

/*static*/ double Aligner::default_delta(bool is_insertion)
{
    return -4;
}

/*static*/ double Aligner::max(const double aa, const double bb, const double cc)
{
    return std::max(aa, std::max(bb, cc));
}

}
