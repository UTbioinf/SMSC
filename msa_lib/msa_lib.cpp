#include <algorithm>
#include "msa_lib.h"

namespace msa_lib
{

MSA::MSA(): 
        epsi(2), scale(2), Con(3), proj(2), 
        dd(2), hh(2), vv(2), 
        Tpointer(3),
        face(3), costs(2),
        avail_vertex(NULL), avail_edge(NULL), avail_coordinate(NULL)
{
    init();
}

bool MSA::is_match_or_mismatch(size_t loc_s1, size_t loc_s2, const int s2Len,
        int gap_open, int gap_extend,
        int match_score, int cur_score, const int* score_table)
{
    if(loc_s1 > 0)
    {
        if(loc_s2 > 0)
            return (score_table[(loc_s1 - 1) * s2Len + loc_s2 - 1] + match_score == cur_score);
        else
        {
            cur_score -= match_score + gap_open;
            if(loc_s1 > 1)
                cur_score -= gap_extend * (loc_s1 - 1);
            return (cur_score == 0);
        }
    }
    else
    {
        if(loc_s2 > 0)
        {
            cur_score -= match_score + gap_open;
            if(loc_s2 > 1)
                cur_score -= gap_extend * (loc_s2 - 1);
            return (cur_score == 0);
        }
        else
        {
            return (match_score == cur_score);
        }
    }
}

bool MSA::is_insertion(size_t loc_s1, size_t loc_s2, const int s2Len,
        int gap_open, int gap_extend,
        int cur_score, const int* score_table)
{
    if(loc_s1 > 0)
    {
        return (score_table[ (loc_s1 - 1) * s2Len + loc_s2 ] + gap_open == cur_score);
    }
    else
    {
        cur_score -= gap_open + gap_open;
        if(loc_s2 > 1)
            cur_score -= gap_extend * (loc_s2 - 1);
        return (cur_score == 0);
    }
}

bool MSA::is_deletion(size_t loc_s1, size_t loc_s2, const int s2Len,
        int gap_open, int gap_extend,
        int cur_score, const int* score_table)
{
    if(loc_s2 > 0)
    {
        return (score_table[ loc_s1 * s2Len + loc_s2 - 1 ] + gap_open == cur_score);
    }
    else
    {
        cur_score -= gap_open + gap_open;
        if(loc_s1 > 1)
            cur_score -= gap_extend * (loc_s1 - 1);
        return (cur_score == 0);
    }
}

void MSA::nw_alignment(const std::string& s1, const std::string& s2,
        int match, int mismatch, int gap_open, int gap_extend)
{
    const int s1Len = s1.length();
    const int s2Len = s2.length();
    parasail_matrix_t* mat = parasail_matrix_create("ATGC", match, mismatch);
    parasail_result_t* aln = parasail_nw_table_diag_32(
            s1.c_str(), s1Len,
            s2.c_str(), s2Len,
            -gap_open, -gap_extend,
            mat);
    M.assign(2, std::string());
    // do trace back
    const int* score_table = aln->score_table;
    size_t loc_s1 = s1Len;
    size_t loc_s2 = s2Len;
    int cur_score, match_score;
    while(loc_s1 > 0 && loc_s2 > 0)
    {
        --loc_s1;   --loc_s2;
        cur_score = score_table[ loc_s1 * s2Len + loc_s2 ];
        match_score = (s1.at( loc_s1 ) == s2.at( loc_s2)) ? match : mismatch;
        if(is_match_or_mismatch( loc_s1, loc_s2, s2Len, gap_open, gap_extend, 
                match_score, cur_score, score_table))
        {
            M[0].push_back( s1.at( loc_s1 ) );
            M[1].push_back( s2.at( loc_s2 ) );
            --loc_s1;   --loc_s2;
        }
        else if(is_insertion(loc_s1, loc_s2, s2Len, gap_open, gap_extend, cur_score, score_table))
        {
            M[0].push_back( s1.at( loc_s1 ) );
            M[1].push_back( '-' );
            --loc_s1;
        }
        else if(gap_open >= gap_extend || is_deletion(loc_s1, loc_s2, s2Len, gap_open, gap_extend,
                cur_score, score_table))
        {
            assert((loc_s2 > 0 && score_table[ loc_s1 * s2Len + loc_s2 - 1 ] + gap_open == cur_score) ||
                (loc_s2 == 0 && gap_open == cur_score));
            M[0].push_back( '-' );
            M[1].push_back( s2.at( loc_s2 ) );
            --loc_s2;            
        }
        else
        {
            size_t tmp_loc_ins = loc_s1 + 1, tmp_loc_del = loc_s2 + 1;
            size_t tmp_mat_loc_ins = (loc_s1 - 1) * s2Len + loc_s2,
                   tmp_mat_loc_del = loc_s1 * s2Len + loc_s2 - 1;
            int ins_score = cur_score, del_score = cur_score;
            int tmp_cnt_ins = 0, tmp_cnt_del = 0;
            bool found_it = false;
            while(tmp_loc_ins > 0 || tmp_loc_del > 0)
            {
                if(tmp_loc_ins > 0)
                {
                    ++tmp_cnt_ins;
                    if(tmp_loc_ins == 1)
                    {
                        if(is_insertion(0, loc_s2, s2Len, gap_open, gap_extend, ins_score, score_table))
                        //if(gap_open == ins_score)
                        {
                            //add_alignment(cur_state, 'I', tmp_cnt_ins);
                            //loc_s1 = tmp_loc_ins - 2;
                            M[1].append( tmp_cnt_ins, '-' );
                            while(loc_s1 > tmp_loc_ins - 2)
                                M[0].push_back( s1.at( loc_s1-- ) );
                            found_it = true;
                            break;
                        }
                        --tmp_loc_ins;
                    }
                    else
                    {
                        if(score_table[ tmp_mat_loc_ins ] + gap_open == ins_score)
                        {
                            //add_alignment(cur_state, 'I', tmp_cnt_ins);
                            //loc_s1 = tmp_loc_ins - 2;
                            M[1].append( tmp_cnt_ins, '-' );
                            while(loc_s1 > tmp_loc_ins - 2)
                                M[0].push_back( s1.at(loc_s1--) );
                            found_it = true;
                            break;
                        }
                        --tmp_loc_ins;
                        tmp_mat_loc_ins -= s2Len;
                        ins_score -= gap_extend;
                    }
                }
                if(tmp_loc_del > 0)
                {
                    ++tmp_cnt_del;
                    if(tmp_loc_del == 1)
                    {
                        if(is_deletion(loc_s1, 0, s2Len, gap_open, gap_extend, del_score, score_table))
                        //if(gap_open == del_score)
                        {
                            //add_alignment(cur_state, 'D', tmp_cnt_del);
                            //loc_s2 = tmp_loc_del - 2;
                            M[0].append( tmp_cnt_del, '-' );
                            while(loc_s2 > tmp_loc_del - 2)
                                M[1].push_back( s2.at(loc_s2--) );
                            found_it = true;
                            break;
                        }
                        --tmp_loc_del;
                    }
                    else
                    {
                        if(score_table[ tmp_mat_loc_del ] + gap_open == del_score)
                        {
                            //add_alignment(cur_state, 'D', tmp_cnt_del);
                            //loc_s2 = tmp_loc_del - 2;
                            M[0].append( tmp_cnt_del, '-' );
                            while(loc_s2 > tmp_loc_del - 2)
                                M[1].push_back( s2.at(loc_s2--) );
                            found_it = true;
                            break;
                        }
                        --tmp_loc_del;
                        --tmp_mat_loc_del;
                        del_score -= gap_extend;

                    }
                
                }
            }
            if(not found_it)
            {
                std::cerr << '[' <<  __FILE__ << ' ' << __LINE__ << "] Ooops! Something wrong!!!" << std::endl;
                std::cerr << "    loc_s1 = " << loc_s1 << std::endl;
                std::cerr << "    loc_s2 = " << loc_s2 << std::endl;
                std::cerr << "    table[] = " << score_table[ loc_s1 * s2Len + loc_s2 ] << std::endl;
                if(loc_s1 > 0)  std::cerr << "    table[loc_s1 - 1, ] = " << score_table[ (loc_s1 - 1) * s2Len + loc_s2] << std::endl;
                if(loc_s2 > 0)  std::cerr << "    table[, loc_s2 - 1] = " << score_table[ loc_s1 * s2Len + loc_s2 - 1] << std::endl;
                if(loc_s1 > 0 && loc_s2 > 0)    std::cerr << "    table[loc_s1 - 1, loc_s2 - 1] =" << score_table[ (loc_s1 - 1) * s2Len + loc_s2 - 1]<< std::endl;
                exit(1);
            }            
        }
        ++loc_s1;   ++loc_s2;
    }
    while(loc_s1 >= 1)
    {
        M[0].push_back( s1.at(--loc_s1) );
        M[1].push_back( '-' );
    }
    while(loc_s2 > 1)
    {
        M[0].push_back( '-' );
        M[1].push_back( s2.at(--loc_s2) );
    }

    parasail_result_free(aln);
    parasail_matrix_free(mat);
    std::reverse(M[0].begin(), M[0].end());
    std::reverse(M[1].begin(), M[1].end());
}

void MSA::init()
{
    S.clear();
    max_len = 0;
    S.push_back("");
}

void MSA::add_seq(const std::string& str)
{
    S.push_back("-"+ str);
    max_len = std::max( max_len, str.length() + 1 );
}

void MSA::set_NUC_matrix(long long match, long long mismatch, long long gap)
{
    GG = G = gap;
    D['A']['T'] = D['A']['G'] = D['A']['C'] = mismatch;
    D['T']['A'] = D['T']['G'] = D['T']['C'] = mismatch;
    D['G']['A'] = D['G']['T'] = D['G']['C'] = mismatch;
    D['C']['A'] = D['C']['T'] = D['C']['G'] = mismatch;
    D['A']['-'] = D['T']['-'] = D['G']['-'] = D['C']['-'] = gap;
    D['-']['A'] = D['-']['T'] = D['-']['G'] = D['-']['C'] = gap;
    D['A']['A'] = D['T']['T'] = D['G']['G'] = D['C']['C'] = match;
    D['-']['-'] = 0;

    T[0][0][0][0] = 0; /* -- : -- */
    T[0][0][0][1] = G; /* -- : -x */
    T[0][0][1][0] = G; /* -x : -- */
    T[0][0][1][1] = 0; /* -x : -x */
    T[0][1][0][0] = 0; /* -- : x- */
    T[0][1][0][1] = 0; /* -- : xx */
    T[0][1][1][0] = G; /* -x : x- */
    T[0][1][1][1] = 0; /* -x : xx */
    T[1][0][0][0] = 0; /* x- : -- */
    T[1][0][0][1] = G; /* x- : -x */
    T[1][0][1][0] = 0; /* xx : -- */
    T[1][0][1][1] = 0; /* xx : -x */
    T[1][1][0][0] = 0; /* x- : x- */
    T[1][1][0][1] = G; /* x- : xx */
    T[1][1][1][0] = G; /* xx : x- */
    T[1][1][1][1] = 0; /* xx : xx */
    T[2][0][2][0] = 0;
    T[2][0][2][1] = 0;
    T[2][1][2][0] = 0;
    T[2][1][2][1] = 0;
    T[0][2][0][2] = 0;
    T[0][2][1][2] = 0;
    T[1][2][0][2] = 0;
    T[1][2][1][2] = 0;
    T[2][2][2][2] = 0;
}

void MSA::run_msa(bool optimal)
{
#ifdef DEBUG_MSA
    std::cerr << "BIG = " << BIG << ", max_len = " << max_len << std::endl;
#endif
    delta = -1;

    dd.set_dimensions( max_len, max_len << 1 );
    hh.set_dimensions( max_len, max_len << 1 );
    vv.set_dimensions( max_len, max_len << 1 );
    size_t K = S.size();
    epsi.set_dimensions(K, K);
    scale.set_dimensions(K, K);
    Con.set_dimensions(K, K, max_len + 1);
    proj.set_dimensions(K, K);
    Tpointer.set_dimensions( K, K, 3 );

    face.set_dimensions( K, K, max_len );
    costs.set_dimensions( K, K );
    msa_A.assign(K + 1, NULL);
    M.assign(K - 1, std::string());

    if(K <= 3) // run parasail
    {
        std::cerr << "Too few strings, call nw instead!" << std::endl;
        exit(1);
    }
    // Calculate pairwise alignments to set epsilons
    primer();
    bias();
    ecalc( 2 * max_len - 2, optimal );

    if(optimal)
    {
        faces();
        if(delta < 0)
        {
            delta = 0;
            for(size_t i = 1; i < K; ++i)
                for(size_t j = i + 1; j < K; ++j)
                    delta += (scale.at(i, j) * epsi.at(i, j));
        }
        Upper = Lower + delta;
        for(size_t i = 1; i < K; ++i)
            for(size_t j = i + 1; j < K; ++j)
            {
                proj.at(i, j) = 0;
                scale.at(j, i) = scale.at(i, j); // make symmetric for msa()
            }
        msa_result = msa();
    #ifdef DEBUG_MSA
        display(msa_result);
    #endif
        compute_opt_aln( msa_result );
        free_msa();
    }
}

void MSA::bias()
{
    MSA_Bias a( *this );
    a.bias();
}

void MSA::ecalc(long long len, bool optimal)
{
    MSA_Ecalc a( *this );
    a.ecalc( len, optimal );
}

void MSA::primer()
{
    size_t K = S.size() - 1;
    for(size_t I = 1; I <= K; ++I)
        for(size_t i = 0; i < S[I].length(); ++i)
            Con.at(I, I, i) = i;
    for(size_t I = 1; I < K; ++I)
    {
        const std::string& A = S[I];
        size_t n = A.length() - 1;
        for(size_t J = I + 1; J <= K; ++J)
        {
            const std::string& B = S[J];
            size_t m = B.length() - 1;

            // compute distance from <0, 0> to <i, j>
            dd.at(0, 0) = 0;
            hh.at(0, 0) = vv.at(0, 0) = GG;
            Con.at(I, J, 0) = Con.at(J, I, 0) = 0;
            for(size_t j = 1; j <= m; ++j)
            {
                vv.at(0, j) = dd.at(0, j) = BIG;
                hh.at(0, j) = hh.at(0, j - 1) + D[DASH][B.at(j)];
                Con.at(J, I, j) = -1;
            }
            for(size_t i = 1; i <=n; ++i)
            {
                hh.at(i, 0) = dd.at(i, 0) = BIG;
                vv.at(i, 0) = vv.at(i-1, 0) + D[A.at(i)][DASH];
                Con.at(I, J, i) = -1;
            }
            for(size_t i = 1; i <= n; ++i)
            {
                long long Gi = i == n ? GG : G;
                for(size_t j = 1; j <=m; ++j)
                {
                    long long Gj = j == m ? GG : G;
                    dd.at(i, j) = min3( dd.at(i-1, j-1), hh.at(i-1, j-1), vv.at(i-1, j-1)) + D[A.at(i)][B.at(j)];
                    hh.at(i, j) = min3( dd.at(i, j-1) + Gi, hh.at(i, j-1), vv.at(i, j-1) + Gi) + D[DASH][B.at(j)];
                    vv.at(i, j) = min3( dd.at(i-1, j) + Gj, hh.at(i-1, j) + Gj, vv.at(i-1, j)) + D[A.at(i)][DASH];
                }
                
            }

            costs.at(I, J) = min3( dd.at(n, m), hh.at(n, m), vv.at(n, m) );
            scale.at(J, I) = convert(I, J, n, m);
        }
    }
}

long long MSA::min3(long long a, long long b, long long c)
{
    if(a < b)
    {
        if(a < c)   return a;
        else    return c;
    }
    else
    {
        if(b < c)   return b;
        else    return c;
    }
}

long long MSA::convert(long long I, long long J, long long n, long long m)
{
    long long i, j, V, H, M;
    long long dir = DIAG;
    long long match = 0;

    for(i = n, j = m; i || j;)
    {
        V = vv.at(i, j) - (dir == VERT ? (j == m ? GG:G) : 0);
        H = hh.at(i, j) - (dir == HORZ ? (i == n ? GG:G) : 0);
        M = min3(V, H, dd.at(i, j));
        if(!j || M == V)
        {
            dir = VERT; --i;
        }
        else if(!i || M == H)
        {
            dir = HORZ; --j;
        }
        else
        {
            dir = DIAG;
            match += S[I].at(i) == S[J].at(j);
            Con.at(I, J, i) = j;
            Con.at(J, I, j) = i;
            --i;    --j;
        }
    }
    return ((long long) (0.5 + 1000.0 * (n + m - match - match) / (n + m)));
}

void MSA::faces()
{
    long long i, j, Gi, Gj;
    std::vector<long long> d_(max_len); // reverse diagonal distance
    std::vector<long long> h_(max_len); // reverse horizontal distance
    std::vector<long long> v_(max_len); // reverse vertical distance
    std::vector<long long> col(max_len);
    long long n, m, U, I, J, w, h_l, h_r, d_l, d_r, v_l, v_r;
    ROW* row;

    size_t K = S.size();
    for(Lower = 0, I = 1; I < K - 1; ++I)
    {
        const std::string& A = S[I];
        n = A.length() - 1;
        for(J = I + 1; J < K; ++J)
        {
            const std::string& B = S[J];
            m = B.length() - 1;

            // compute distance from <0, 0> to <i, j>
            dd.at(0, 0) = 0; hh.at(0, 0) = vv.at(0, 0) = GG;
            for(j = 1; j <= m; ++j)
            {
                vv.at(0, j) = dd.at(0, j) = BIG;
                hh.at(0, j) = hh.at(0, j - 1) + D[ DASH ][ B.at(j) ];
            }
            for(i = 1; i <= n; ++i)
            {
                hh.at(i, 0) = dd.at(i, 0) = BIG;
                vv.at(i, 0) = vv.at(i-1, 0) + D[ A.at(i) ][ DASH ];
            }
            for(i = 1; i <= n; ++i)
            {
                Gi = i == n ? GG : G;
                for(j = 1; j <= m; ++j)
                {
                    Gj = j == m ? GG : G;
                    dd.at(i, j) = min3(dd.at(i-1, j-1), hh.at(i-1, j-1), vv.at(i-1, j-1)) + D[A.at(i)][B.at(j)];
                    hh.at(i, j) = min3(dd.at(i, j-1) + Gi, hh.at(i, j-1), vv.at(i,j-1) + Gi) + D[DASH][B.at(j)];
                    vv.at(i, j) = min3(dd.at(i-1, j) + Gj, hh.at(i-1, j) + Gj, vv.at(i-1, j)) + D[A.at(i)][DASH];
                }
            }
            U = (costs.at(I, J) = min3(dd.at(n, m), hh.at(n, m), vv.at(n, m))) + epsi.at(I, J);
            Lower += scale.at(I, J) * costs.at(I, J);

            // compute distance from <n, m> to <i, j>
            d_.at(m) = 0; h_.at(m) = v_.at(m) = GG;
            for( j = m - 1; j >= 0; --j)
                v_.at(j) = (d_.at(j) = h_.at(j) = h_.at(j+1) + D[DASH][B.at(j+1)]) + G;
            for( j = w = 0; j < m; ++j)
                if(min3(hh.at(n, j) - GG, dd.at(n, j), vv.at(n, j)) + h_.at(j) <= U)
                    col.at(w++) = j;
            col.at(w++) = m;

            row = &face.at(I, J, n);
            row->width = w;
            row->column.assign(col.begin(), col.begin() + w);
            for(i = n - 1; i >= 0; --i)
            {
                Gi = i == 0 ? GG : G;
                h_r = h_.at(m);    d_r = d_.at(m);    v_r = v_.at(m);
                h_.at(m) = (d_.at(m) = v_.at(m) = v_r + D[A.at(i+1)][DASH]) + G;
                for(j = m - 1; j >= 0; --j)
                {
                    Gj = j == 0 ? GG : G;
                    h_l = h_.at(j);    d_l = d_.at(j);    v_l = v_.at(j);
                    d_.at(j) = min3(d_r, h_r, v_r) + D[A.at(i+1)][B.at(j+1)];
                    h_.at(j) = min3(d_.at(j+1) + Gi, h_.at(j+1), v_.at(j+1) + Gi) + D[DASH][B.at(j+1)];
                    v_.at(j) = min3(d_l + Gj, h_l + Gj, v_l) + D[A.at(i+1)][DASH];
                    h_r = h_l;  d_r = d_l;  v_r = v_l;
                }
                for(j = w = 0; j <= m; ++j)
                {
                    Gj = j == 0 || j == m ? GG : G;
                    if(min3(
                            hh.at(i, j) + min3(h_.at(j) - Gi, d_.at(j), v_.at(j)),
                            dd.at(i, j) + min3(h_.at(j)     , d_.at(j), v_.at(j)),
                            vv.at(i, j) + min3(h_.at(j)     , d_.at(j), v_.at(j) - Gj)) <= U)
                        col.at( w++ ) = j;
                }
                row = &face.at(I, J, i);
                row->width = w;
                row->column.assign(col.begin(), col.begin() + w);
            }
        }
    }
}

void MSA::INSERT(EDGE* &e, HEAP &h)
{
    EDGE* &b = h.bucket.at( e->dist );
    if(b != NULL)
        b->heap_pred = e;
    e->heap_succ = b;
    e->heap_pred = NULL;
    b = e;
}

void MSA::DELETE(EDGE* &e, HEAP &h)
{
    EDGE* &b = h.bucket.at( e->dist );
    if(e->heap_pred != NULL)
        e->heap_pred->heap_succ = e->heap_succ;
    else
        b = e->heap_succ;
    if(e->heap_succ != NULL)
        e->heap_succ->heap_pred = e->heap_pred;
}

EDGE* MSA::msa()
{
    size_t K = S.size();
    std::vector<long long> p(K), q(K), r(K);
    long long d;
    VERTEX *v, *t;
    VERTEX *s;
    EDGE *e, *f;
    std::string C(K, 0);
    long long I, J;
    std::vector<long long> delta0(K), delta1(K);
    long long difference;
    std::vector<long long> ends(K-1);

    // compute shortest paths to vertices in intersected region of lattice
    s = source();
    t = sink();
    heap(Upper);
    presource = create_vertex(NULL);
    e = create_edge( presource, s );
    e->dist = 0;
    e->refer++; // make sure edge does not get freed
    e->backtrack = NULL;
    INSERT(e, h);

    while ( (e= extract()) != NULL && (v=e->head) != t)
    {
        if(e->dist <= Upper)
        {
            // put coordinates of tail into p array
            // and coordinates of head of edge into q
            coord(e->tail, p);
            safe_coord(e->head, q);
            // next loop is from cost function
            // difference between p and q
            for(I = 1; I < K; ++I)
                delta0.at(I) = q.at(I) - p.at(I);
            for(I = 2; I < K; ++I)
                for(J = 1; J < I; ++J)
                {
                    Tpointer.at(I, J, 0) = T[ delta0.at(I) ][ delta0.at(J) ][ 0 ];
                    Tpointer.at(I, J, 1) = T[ delta0.at(I) ][ delta0.at(J) ][ 1 ];
                    Tpointer.at(I, J, 2) = T[ delta0.at(I) ][ delta0.at(J) ][ 2 ];
                }
            if(v->out == NULL) // If first time visiting v
                adjacent(e, q);
            for(f = v->out; f!= NULL; f=f->next)
            {
                difference = f->dist - e->dist;
                if(difference > 0)
                {
                    // get coordinates of next vertex r
                    safe_coord(f->head, r);

                    for(I = 1; I < K; ++I)
                    {
                        C.at(I) = (r.at(I) > q.at(I) ? S.at(I).at(r.at(I)) : DASH);
                        delta1.at(I) = r.at(I) - q.at(I);
                    }
                    d = 0;
                    d += scale.at(1, 2) * ( D[C.at(1)][C.at(2)] + Tpointer.at(2, 1, delta1.at(2))[delta1.at(1)] );
                    for(I = K -1; I >= 3; --I)
                    {
                        if(d >= difference)
                            goto nextedge;
                        for(J = 1; J < I; ++J)
                            d += scale.at(I, J) * (D[C.at(I)][C.at(J)] + Tpointer.at(I, J, delta1.at(I))[delta1.at(J)]);
                    }

                    if(d < difference)
                    {
                        DELETE(f, h);
                        e->refer++;
                        if(f->backtrack != NULL)
                            if( -- f->backtrack->refer == 0)
                                free_edge( f->backtrack );
                        f->dist = d + e->dist;
                        f->backtrack = e;
                        INSERT(f, h);
                    }
                }
                nextedge: continue;
            }
        }
        if(e->refer == 0)
            free_edge( e );
    }
    return e;
}

// create source vertex of lattice
VERTEX* MSA::source()
{
    size_t K = S.size();
    std::vector<long long> p(K, 0), index(max_len);
    long long i;
    COORDINATE *a;
    for(i = S[1].length() - 1; i >=0; --i)
        index.at(i) = i;

    a = msa_A.at(1) = create_coordinate(index, S[1].length(), NULL);
    for(i = 2; i < K; ++i)
        a = a -> coord_vals -> next_coord = create_coordinate(index, intersect(p, i, index), a->coord_vals);
    msa_A.at(1)->refer++;
    return (VERTEX *)( a->coord_vals -> next_coord = (COORDINATE *)create_vertex(a->coord_vals));
}

// create sink vertex of lattice
VERTEX* MSA::sink()
{
    size_t K = S.size();
    long long i;
    std::vector<long long> p(K), index(max_len);
    COORDINATE *a;
    COORDINATE_VALUES *f;
    for(i = 1; i < K; ++i)
        p.at(i) = S.at(i).length() - 1;
    a = msa_A.at(1);
    for(i = 2; i < K; ++i)
    {
        f = a->coord_vals + p.at(i-1) - a->lo;
        a = f->next_coord = create_coordinate(index, intersect(p, i, index), f);
    }
    f = a->coord_vals + p.at(K-1) - a->lo;
    return (VERTEX*)(f->next_coord = (COORDINATE *)create_vertex(f));
}

// intersect regions on rows of faces
// Finds the possible values for the "seqnum"th coordinate of point,
// given point[1], point[2], ..., point[seqnum - 1]. The cooresponding
// bound values are copied to bound[].
long long MSA::intersect(std::vector<long long>& point, long long seqnum, std::vector<long long>& possible_values)
{
    long long i, j, k;
    long long J, m, n;
    ROW* r;

    // retrieve values that are consistent with the pairwise alignment of sequence 1 and seqnum
    r = &face.at(1, seqnum, point.at(1));
    m = r->width;
    for(i = 0; i < m; ++i)
        possible_values.at( i ) = r->column.at(i);
    for(J = 2; J < seqnum; ++J)
    {
        // Get values that are consistent with the pairwise
        // alignment of sequences J and seqnum
        r = &face.at(J, seqnum, point.at(J));
        n = r->width;
        // compute intersection of possible values array and
        // c array for sequence j, keeping result in possible
        // values
        for(i = j = k = 0; i < m && j < n; )
        {
            if(possible_values.at(i) < r->column.at(j))
                ++i;
            else if(possible_values.at(i) > r->column.at(j))
                ++j;
            else
            {// possible_values[i] == c[j]
                possible_values.at(k++) = possible_values.at(i++);
                ++j;
            }
        }
        if((m = k) <= 0)    break;
    }
    return m;
}

// generate adjacent vertices within region
void MSA::adjacent(EDGE* e, std::vector<long long>& q)
{
    size_t K = S.size();
    long long I, i;
    std::vector<long long> index( max_len );
    long long n;
    COORDINATE *a;
    COORDINATE_VALUES *f;
    std::vector<long long> P(K), B(K);

    for(i = K - 1, f = e->head->prev_coord_val; i >= 1; --i, f = f->curr_coord->prev_coord_val)
    {
        msa_A.at(i) = f->curr_coord;
        B.at(i) = (P.at(i) = q.at(i)) + 1;
    }

    for(I = K - 1; I < K; )
    {
        for(; I >= 1; --I)
        {
            if(P.at(I) < B.at(I))
                break;
        }
        if(I < 1)
            return;
        // first test is to see whether i is in the proper range,
        // second test is to see whether i is a valid coordinate
        if( (i = P.at(I) = B.at(I)) >= (a = msa_A.at(I))->lo && i <= a->hi &&
                (f = a->coord_vals + i - a->lo)->curr_coord != NULL)
        {// if child in trie does not exist
            if(f->next_coord == NULL)
                if(I < K - 1)
                    if( (n = intersect(P, I+1, index)) > 0)
                        f->next_coord = create_coordinate(index, n, f);
                    else
                    {
                        f->curr_coord = NULL;
                        if(free_coordinate(a))
                            --I;
                        continue;
                    }
                else    f->next_coord = (COORDINATE *)create_vertex(f);
            for(a = msa_A.at(++I) = f->next_coord; I < K; a = msa_A.at(++I) = f->next_coord)
                if((i = P.at(I) = B.at(I) - 1) >= a->lo && i <= a->hi &&
                        (f = a->coord_vals + i - a->lo) -> curr_coord != NULL)
                {
                    if(f->next_coord == NULL)
                        if( I < K-1)
                            if( (n = intersect(P, I+1, index)) > 0 )
                                f->next_coord = create_coordinate(index, n, f);
                            else
                            {
                                f->curr_coord = NULL;
                                if( free_coordinate(a))
                                    --I;
                                break;
                            }
                        else    f->next_coord = (COORDINATE *)create_vertex(f);
                }
                else    break;
        }
        if(I == K)
        {
            create_edge(e->head, (VERTEX *)a);
            I = K - 1;
        }
    }
}

// compute lattice coordinate of vertex
void MSA::safe_coord(VERTEX* v, std::vector<long long>& p)
{
    long long i;
    COORDINATE_VALUES *a;
    COORDINATE *b;
    for(i = S.size() - 1, a = v->prev_coord_val; i >= 1; --i, a = b->prev_coord_val)
    {
        b = a->curr_coord;
        p.at(i) = a - b->coord_vals + b->lo;
    }
}

// compute lattice coordinate of vertex
void MSA::coord(VERTEX* v, std::vector<long long>& p)
{
    long long i;
    COORDINATE_VALUES *a;
    COORDINATE *b;
    if( v!=presource)
    {
        for(i = S.size() - 1, a = v->prev_coord_val; i >= 1; --i, a=b->prev_coord_val)
        {
            b = a->curr_coord;
            p.at(i) = a - b->coord_vals + b->lo;
        }
    }
    else
    {
        for(i = S.size() - 1; i >= 1; --i)
            p.at(i) = -1;
    }
}

// compute projected cost of edge <q, r> preceded by <p, q>
void MSA::project(std::vector<long long>& p, std::vector<long long>& q, std::vector<long long>& r)
{
    size_t K = S.size();
    std::string C(K, 0);
    long long I, J;
    std::vector<long long> t0(K), t1(K);
    for(I = 1; I < K; ++I)
    {
        C.at(I) = (r.at(I) > q.at(I) ? S.at(I).at( r.at(I) ) : DASH);
        t0.at(I) = q.at(I) - p.at(I);
        t1.at(I) = r.at(I) - q.at(I);
    }
    for(I = 1; I < K; ++I)
        for(J = I + 1; J < K; ++J)
            proj.at(I, J) += D[C.at(I)][C.at(J)] + T[ t0.at(I) ][ t0.at(J) ][ t1.at(I) ][ t1.at(J) ];
}

void MSA::compute_opt_aln(EDGE* e)
{
    size_t K = S.size();
    std::vector<long long> p(K), q(K), r(K);
    long long d;
    EDGE* f;

    if(e == NULL || (d = e->dist) > Upper)
    {
        std::cerr << "[" << __FILE__ << ' ' << __LINE__ << "]: Multiple alignment within bound does not exist" << std::endl;
        exit(1);
    }

    for(e->next = NULL; e ->tail != presource; e = f)
    {
        f = e->backtrack;
        coord( f->tail, p );
        safe_coord(f->head, q);
        safe_coord(e->head, r);
        f->next = e;
        project(p, q, r);
    }
    for(e = e->next; e != NULL; e = e->next)
    {
        coord(e->tail, p);
        safe_coord(e->head, q);
        for(size_t k = 1; k < K; ++k)
            M.at(k - 1).push_back( S.at(k).at( q.at(k)*(q.at(k) - p.at(k)) ) );
    }
}

#ifdef DEBUG_MSA
void MSA::display(EDGE* e)
{
    C = 0;
    size_t K = S.size();
    std::vector<long long> p(K), q(K), r(K);
    long long d;
    EDGE* f;
    // shortest sink to source path in lattice within bound
    if(e == NULL || (d=e->dist) > Upper)
    {
        std::cerr << "[" << __FILE__ << ' ' << __LINE__ << "]: Multiple alignment within bound does not exist" << std::endl;
        exit(1);
    }
    // recover shortest path to source tracing backward from sink
    for(e->next = NULL; e ->tail != presource; e=f)
    {
        f = e->backtrack;
        coord( f->tail, p );
        safe_coord(f->head, q);
        safe_coord(e->head, r);
        f->next = e;
        project(p, q, r);
    }
    // display alignment corresponding to path
    std::cout << "\n                  ***  Optimal Multiple Alignment  ***\n\n";
    M.resize( S.size(), std::string(LINE, 0) );
    for(e = e->next; e!=NULL; e=e->next)
    {
        coord(e->tail, p);
        safe_coord(e->head, q);
        column(p, q);
    }
    output();
    // output statistics
    /// std::cout << "Alignment cost: " << d << "    Lower bound: " << Lower << std::endl;
    /// std::cout << "Delta:          " << d-Lower << "    Max. Delta: " << delta << std::endl << std::endl;
    /// std::cout << "Sequences  Proj. Cost  Pair. Cost  Epsilon";
    /// std::cout << "  Max. Epsi.  Weight  Weight*Cost" << std::endl;
    /// for(I = 1; I < K; ++I)
    ///     for(J = I + 1; J < K; ++J)
    ///         std::cout << I << '\t' << J << '\t' << proj.at(I, J) << '\t' << costs.at(I, J) << '\t'
    ///                   << proj.at(I, J) - costs.at(I, J) << '\t' << epsi.at(I, J) << '\t' 
    ///                   << scale.at(I, J) << '\t' << scale.at(I, J) * proj.at(I, J) << std::endl;
}

void MSA::column(std::vector<long long>& p, std::vector<long long>& q)
{
    long long k;
    for(k = 1; k < S.size(); ++k)
        M.at(k).at(C) = S.at(k).at( q.at(k) * (q.at(k) - p.at(k)));
    if(++C >= LINE)
        output();
        
}

void MSA::output()
{
    long long k, c;
    if(C == 0)  return;
    for(k = 1; k < S.size(); ++k)
    {
        for(c = 0; c < C; ++c)
            std::cout << M.at(k).at(c);
        std::cout << std::endl;
    }
    std::cout << std::endl;
    C = 0;
}
#endif

void MSA::heap(long long max)
{
    h.bucket.assign( max + 2, NULL );
    h.min = 0;
    h.max = max;
}

// extract minimum distance edge from heap
EDGE* MSA::extract()
{
    EDGE *e;
    for(size_t i = h.min; i <= h.max; ++i)
        if(h.bucket.at(i) != NULL)
        {
            if( (h.min = i) > h.max )
                return NULL;
            e = h.bucket.at(i);
            break;
        }
    DELETE(e, h);

    // remove e from the list of non-extracted in-edges incident to the head of e
    if(e->nonextracted_prev != NULL)
        e->nonextracted_prev->nonextracted_next = e->nonextracted_next;
    else
        e->head->nonextracted_inedges = e->nonextracted_next;
    if(e->nonextracted_next != NULL)
        e->nonextracted_next->nonextracted_prev = e->nonextracted_prev;

    e->nonextracted_prev = NULL;
    e->nonextracted_next = NULL;
    return e;
}

VERTEX* MSA::create_vertex(COORDINATE_VALUES* prev_coord_val)
{
    VERTEX *v;
    if(avail_vertex != NULL)
    {
        v = avail_vertex;
        avail_vertex = (VERTEX *)v -> out;
    }
    else
    {
        v = new VERTEX();
        vertex_pointers.push_back( v );
    }
    v->prev_coord_val = prev_coord_val;
    v->out = NULL;
    v->nonextracted_inedges = NULL;
    return v;
}

void MSA::free_vertex(VERTEX* v)
{
    COORDINATE *a;
    EDGE *e, *e2;
    a = v->prev_coord_val->curr_coord;
    v->prev_coord_val->curr_coord = NULL;
    free_coordinate(a);
    // remove all the remaining edges coming into v from the heap
    e = v->nonextracted_inedges;
    while(e != NULL)
    {
        DELETE(e, h);
        e2 = e->nonextracted_next;
        e->refer = -1; // kludge telling free_edge() not to call free_vertex()

        free_edge(e);
        e = e2;
    }
    v->out = (EDGE *) avail_vertex;
    avail_vertex = v;
}

EDGE* MSA::create_edge(VERTEX* v, VERTEX* w)
{
    EDGE *e;
    if(avail_edge != NULL)
    {
        e = avail_edge;
        avail_edge = e->next;
    }
    else
    {
        e = new EDGE();
        edge_pointers.push_back( e );
    }
    // set endpoints of e
    e->tail = v;
    e->head = w;
    // insert e into v's list of outgoing edges
    e->prev = NULL;
    e->next = v->out;
    e->backtrack = NULL;
    e->refer = 0;
    v->out = e;
    if(e->next != NULL)
        e->next->prev = e;
    e->heap_succ = e->heap_pred = e->nonextracted_next = e->nonextracted_prev = e;
    e->dist = Upper + 1; // set e's distance to upper bound for lack of more info
    // insert e into w's list of nonextracted incoming edges
    if(e->head->nonextracted_inedges != NULL)
        e->head->nonextracted_inedges->nonextracted_prev = e;
    e->nonextracted_next = e->head->nonextracted_inedges;
    e->nonextracted_prev = NULL;
    e->head->nonextracted_inedges = e;
    return e;
}

void MSA::free_edge(EDGE* e)
{
    if(e->heap_succ != NULL && e->heap_pred != NULL)
        DELETE(e, h);
    if(e->nonextracted_prev != NULL && e->nonextracted_next != NULL)
    {
        if(e->nonextracted_prev != NULL)
            e->nonextracted_prev->nonextracted_next = e->nonextracted_next;
        else
            e->head->nonextracted_inedges = e->nonextracted_next;
        if(e->nonextracted_next != NULL)
            e->nonextracted_next->nonextracted_prev = e->nonextracted_prev;
    }

    if(e->prev != NULL) e->prev->next = e->next;
    else e->tail->out = e->next;
    if(e->next != NULL) e->next->prev = e->prev;
    if((e->backtrack) && (-- e->backtrack->refer == 0))
        free_edge(e->backtrack);
    if((e->head->out == NULL) && (e->refer != -1))
        free_vertex(e->head);
    e->next = avail_edge;
    avail_edge = e;
}

COORDINATE* MSA::create_coordinate(std::vector<long long>& index, long long n, COORDINATE_VALUES* prev_coord_val)
{
    COORDINATE_VALUES *p;
    COORDINATE *a;
    // look for available COORDINATE on free list
    if( avail_coordinate != NULL)
    {
        a = avail_coordinate;
        avail_coordinate = a->next_on_free_list;
    }
    else
    {
        a = new COORDINATE();
        coordinate_pointers.push_back( a );
        a->next_on_free_list = NULL;
    }
    a->lo = index.at(0);
    a->hi = index.at(n-1);
    a->prev_coord_val = prev_coord_val;
    a->refer = n;
    a->coord_vals = new COORDINATE_VALUES[a->hi - a->lo + 1];
    for(p = a->coord_vals + a->hi - a->lo; p >= a->coord_vals; p--)
        p->next_coord = p->curr_coord = NULL;
    // link each new COORDINATE_VALUE it its parent COORDINATE and adjust
    // bounds; note that this loop works only on the possible values which
    // are stored in the array index
    for(n--; n>=0; n--)
        (p = a->coord_vals + index.at(n) - a->lo)->curr_coord = a;
    return a;
}

bool MSA::free_coordinate(COORDINATE *a)
{
    COORDINATE* b;
    if( -- a->refer <= 0)
    {
        b = a->prev_coord_val->curr_coord;
        a->prev_coord_val->curr_coord = NULL;
        delete[] (a->coord_vals);
        a->coord_vals = NULL;
        a->next_on_free_list = avail_coordinate;
        avail_coordinate = a;
        if(b != NULL)
            free_coordinate(b);
        return true;
    }
    else return false;
}

void MSA::free_msa()
{
    // free memory for the data structure
    // free memory for vertices
    for(std::vector<VERTEX *>::iterator it = vertex_pointers.begin();
            it != vertex_pointers.end(); ++it)
        delete *it;
    avail_vertex = NULL;
    vertex_pointers.clear();

    // free memory for edges
    for(std::vector<EDGE *>::iterator it = edge_pointers.begin();
            it != edge_pointers.end(); ++it)
        delete *it;
    avail_edge = NULL;
    edge_pointers.clear();

    // free memory for coordinates
    for(std::vector<COORDINATE *>::iterator it = coordinate_pointers.begin();
            it != coordinate_pointers.end(); ++it)
    {
        if((*it)->coord_vals != NULL)
            delete[] ((*it)->coord_vals);
        delete *it;
    }
    avail_coordinate = NULL;
    coordinate_pointers.clear();
}

std::vector<std::string>& MSA::get_alignment()
{
    return M;
}

}
