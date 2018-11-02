#include "msa_ecalc.h"

namespace msa_lib
{

MSA_Ecalc::MSA_Ecalc(MSA& m): msa(m),
        AL(2), NAL(2), OAL(2)
{
    K = msa.S.size();
}

void MSA_Ecalc::ecalc(long long len, bool optimal)
{
    SS.assign(K, NULL);
    AL.set_dimensions(K, msa.max_len <<1);
    NAL.set_dimensions(K, msa.max_len << 1);
    OAL.set_dimensions(K, msa.max_len << 1);
    Ord.assign(K, 0);

    long long I, J, i, j, c , pos, p1, oldp1, test;

    for(I = 1; I < K; ++I)
        for(J = 1; J < K; ++J)
            msa.Con.at(I, J, msa.S.at(I).length()) = msa.S.at(J).length();
    for(I = 1; I < K; ++I)
        AL.at(I, 0) = DASH;
    pos = oldp1 = 0;
    // previous parts are tested

    for(p1 = 1; p1 <= msa.S.at(1).length(); ++p1)
    {
        // Search for positions consistent among all pairwise alignments and
        // align segments of all sequences from last consistent position to
        // current one
        for(I = test = 2; test && I < K; ++I)
        {
            test = ((i = msa.Con.at(1, I, p1)) >= 0);
            for(J = I + 1; test && J < K; ++J)
                test = ((j = msa.Con.at(1, J, p1)) >= 0 && msa.Con.at(I, J, i) == j);
        }
        if(test)
        {
            // Pick an order in which to construct a progressive multiple alignment
            order(oldp1, p1);
            // Align all segments in the order chosen
            j = align(oldp1, p1, len);
            // Add the aligned segments to the heuristic alignment
            for(I = 1; I < K; ++I)
            {
                J = Ord.at(I);
                for(i = 1; i <= j; ++i)
                    AL.at(J, pos+i) = msa.S.at(J).at( NAL.at(I, j-i) );
            }
            // record whether a given alignment position is consistent among lal
            // optimal alignments, or has arisen by a progressive alignment

            pos += j + 1;
            for(I = 1; I < K; ++I)
            {
                if(msa.Con.at(1, I, p1) == msa.S.at(I).length())
                    continue;
                AL.at(I, pos) = msa.S.at(I).at( msa.Con.at(1, I, p1) );
            }
            oldp1 = p1;

        }
    }

    --pos;
    if(optimal)
    {
        // for each pair, calculate difference between imposed alignment cost 
        // and optimal alignment cost
        for(I = 1; I < K ; ++I)
            for(J = I + 1; J < K; ++J)
            {
                c = pcost(I, J, pos) - msa.costs.at(I, J);
                msa.epsi.at(I, J) = (c < MINE) ? MINE : (c > MAXE ? MAXE : c);
            }
    }
    else
    {
        // print out the heuristic multiple alignment
        for(I = 1; I < K; ++I)
        {
            for(j = 1; j <= pos; ++j)
                msa.M.at(I - 1).push_back( AL.at(I, j) );
        }

    #ifdef DEBUG_MSA
        std::cout << "\n                 ***  Heuristic Multiple Alignment  ***\n\n";
        for(i = pos; i > 0; i-=LINE)
        {
            for(I = 1; I < K; ++I)
            {
                for(j=1; j <= LINE && j <= i; ++j)
                    std::cout << char(AL.at(I, pos - i + j));
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    #endif
    }

    
}

void MSA_Ecalc::order(long long oldp1, long long p1)
{
    long long I, J, i, j, i1, j1, PI, PJ, t, mind, p, M, end;

    std::vector<long long> test(K, 1);
    std::vector<long long> sum(K);
    NDArray<long long> dis(2);
    dis.set_dimensions(K, K);
    std::string a1, a2;
    a1.reserve( msa.max_len );
    a2.reserve( msa.max_len );

    a1.push_back( DASH );
    a2.push_back( DASH );
    mind = BIG;

    // Calculate all pairwise costs for the segments in question
    for(I = 1; I < K; ++I)
        for(J = I + 1; J < K; ++J)
        {
            PJ = msa.Con.at(1, J, oldp1);
            PI = msa.Con.at(1, I, oldp1);
            end = msa.Con.at(1, I, p1);
            for(p = 0, i = PI + 1; i <= end; ++i)
            {
                j = (t = msa.Con.at(I, J, i-1)) >= 0 ? msa.Con.at(I, J, i) - t : msa.Con.at(I, J, i) > 0;
                for(M = j > 0 ? j : 1; M; --M)
                {
                    if(PI >= msa.S.at(I).length() || PJ >= msa.S.at(J).length())
                    {
                        std::cerr << "[" << __FILE__ << ' ' << __LINE__ << "] [ERROR]: Segment fault" << std::endl;
                        exit(1);
                    }
                    if(++p < a1.length())
                    {
                        if(PI + 1 == msa.S.at(I).length() && M == 1 || PJ + 1 == msa.S.at(J).length() && M <= j)
                            continue;
                        a1[ p ] = M == 1 ? msa.S.at(I).at(++PI) : DASH;
                        a2[ p ] = M <= j ? msa.S.at(J).at(++PJ) : DASH;
                    }
                    else
                    {
                        if(PI + 1 == msa.S.at(I).length() && M == 1 || PJ + 1 == msa.S.at(J).length() && M <= j)
                            continue;
                        a1.push_back( M == 1 ? msa.S.at(I).at(++PI) : DASH );
                        a2.push_back( M <= j ? msa.S.at(J).at(++PJ) : DASH );
                    }
                }
            }
            for(t = 0, i = 1; i < p; ++i)
                t += msa.D[ a1[i] ][ a2[i] ] + msa.T[ a1[i-1]!=DASH ][ a2[i-1]!=DASH ][ a1[i]!=DASH ][ a2[i]!=DASH ];
            if(t < mind)
            {
                mind = t;
                i1 = I;
                j1 = J;
            }
            dis.at(I, J) = dis.at(J, I) = t;
        }
    // Make lowest cost pair the first two segments in the order
    // for(I = 1; I < K; ++I)  test[ I ] = 1; // done by initialization
    Ord.at(1) = i1;    SS.at(1) = &(msa.S.at(i1));  test.at(i1) = 0;
    Ord.at(2) = j1;    SS.at(2) = &(msa.S.at(j1));  test.at(j1) = 0;
    for(I = 1; I < K; ++I)
        sum.at(I) = dis.at(I, i1) + dis.at(I, j1);

    // Fill out the other using average distances
    for(j = 3; j < K; ++j)
    {
        mind = BIG;
        for(I = 1; I < K; ++I)
            if(test.at(I) && sum.at(I) < mind)
                mind = sum.at(i = I);
        Ord.at(j) = i; SS.at(j) = &(msa.S.at(i));   test.at( i ) = 0;
        for(I = 1; I < K; ++I)
            sum.at(I) += dis.at(I, i);
    }
}

// Subroutine to calculate an optimal progressive
// alignment of the segments in the order chosen
long long MSA_Ecalc::align(long long oldp1, long long p1, long long len)
{
    long long I, J, i, low, m;
    for(I = 1; I < K; ++I)
        OAL.at(I, 0) = 0;
    low = msa.Con.at(1, Ord.at(1), oldp1);
#ifdef DEBUG_MSA
    A = SS.at(1)->substr(low);
#else
    A = SS.at(1)->c_str() + low;
#endif
    m = msa.Con.at(1, Ord.at(1), p1) - low - 1;
    for(i = 1; i <=m; ++i)
        NAL.at(1, m-i) = low + i;
    for(I = 2; I < K; ++I)
    {
        for(J = 1; J < I; ++J)
            for(i = 1; i <= m; ++i)
                OAL.at(J, i) = NAL.at(J, m-i);
        low = msa.Con.at(1, Ord.at(I), oldp1);
    #ifdef DEBUG_MSA
        A = SS.at(I)->substr(low);
    #else
        A = SS.at(I)->c_str() + low;
    #endif
        m = newal(I, low, msa.Con.at(1, Ord.at(I), p1) - low - 1, m);
        if(m > len && I < K - 1)
        {
            std::cerr << "[" << __FILE__ << ' ' << __LINE__ << "]: Heuristic alignment segment is longer than " << len << std::endl;
            exit(1);
        }
    }
    return m;
}

// Subroutine to calculate the optimal alignment of a new segment and
// a multiple alignment
long long MSA_Ecalc::newal(long long I, long long low, long long n, long long m)
{
    if(n <=0 || m <= 0) return 0;
    long long i, j, k, x, q, qq, dg, vg, hg, sum, sum2;
    std::vector<long long> d( msa.max_len << 1), v( msa.max_len << 1), h( msa.max_len  << 1);
    std::vector<long long> dn( msa.max_len << 1), vn( msa.max_len << 1 ), hn( msa.max_len << 1 );
    std::vector<long long> sc( K );
    for(sum = sum2 = 0, k = 1; k < I; ++k)
    {
        sum += sc.at(k) = Ord.at(k) < Ord.at(I) ?
                        msa.scale.at(Ord.at(k), Ord.at(I)) : msa.scale.at( Ord.at(I), Ord.at(k) );
        if( OAL.at(k, 1) )
            sum2 += sc.at(k);
    }
    for(j = 0; j <= m; ++j)
        hn.at(j) = dn.at(j) = vn.at(j) = BIG;
    for(i = 0; i <= n; ++i)
    {
        for(j = 0; j <= m; ++j)
        {
            h.at(j) = hn.at(j);   d.at(j) = dn.at(j);   v.at(j) = vn.at(j);
        }
        if( i == 0)
        {
            dn.at(0) = 0;
            vn.at(0) = sum * (low ? msa.G : msa.GG);
            hn.at(0) = sum2 * (low ? msa.G : msa.GG);
        }
        else
        {
        #ifdef DEBUG_MSA
            vn.at(0) = v.at(0) + sum * msa.D[ A.at(i) ][ DASH ];
        #else
            vn.at(0) = v.at(0) + sum * msa.D[ A[i] ][ DASH ];
        #endif
            hn.at(0) = dn.at(0) = BIG;
        }
        msa.vv.at(i, 0) = msa.dd.at(i, 0) = msa.hh.at(i, 0) = 1;

        // Calculate optimal alignment in the case that terminal gap costs
        // are different than internal gap costs
        // ...
        // This part is omitted because end gaps should also be penalized

        // Calculate optimal alignment in the case that terminal gap costs
        // are the same as internal gap costs
        for(j = 1; j <=m ; ++j)
        {
            dg = d.at(j-1);    hg=h.at(j-1);
            for(x = 0, k = 1; k < I; ++k)
            {
                q = OAL.at(k, j - 1);
                qq = OAL.at(k, j);
                dg += sc.at(k) * msa.T[1][q>0][1][qq>0];
                hg += sc.at(k) * msa.T[0][q>0][1][qq>0];
            #ifdef DEBUG_MSA
                x  += sc.at(k) * msa.D[ SS.at(k)->at(qq) ][ A.at(i) ];
            #else
                x  += sc.at(k) * msa.D[ SS.at(k)->at(qq) ][ A[i] ];
            #endif
            }
            dn.at(j) = x + (k = msa.min3(dg, hg, v.at(j-1)));
            msa.dd.at(i, j) = k == dg ? DIAG : (k == hg ? HORZ : VERT);

            dg = d.at(j);  hg = h.at(j);
            for(k = 1; k < I; ++k)
            {
                qq = OAL.at(k, j);
                dg += sc.at(k) * msa.T[1][qq > 0][1][0];
                hg += sc.at(k) * msa.T[0][qq > 0][1][0];
            }
        #ifdef DEBUG_MSA
            vn.at(j) = sum * msa.D[DASH][ A.at(i) ] + (k = msa.min3(dg, hg, v.at(j)));
        #else
            vn.at(j) = sum * msa.D[DASH][ A[i] ] + (k = msa.min3(dg, hg, v.at(j)));
        #endif
            msa.vv.at(i, j) = k == dg ? DIAG : (k == hg ? HORZ : VERT);

            dg = dn.at(j-1);   vg = vn.at(j-1);   hg=hn.at(j-1);
            for(x = 0, k = 1; k < I; ++k)
            {
                q = OAL.at(k, j - 1);
                qq = OAL.at(k, j);
                dg += sc.at(k) * msa.T [1] [q>0] [0] [qq>0];
                vg += sc.at(k) * msa.T [1] [0]   [0] [qq>0];
                hg += sc.at(k) * msa.T [0] [q>0] [0] [qq>0];
                x  += sc.at(k) * msa.D [SS.at(k)->at(qq)] [DASH];
            }
            hn.at(j) = x + (k = msa.min3(dg, hg, vg));
            msa.hh.at(i, j) = k == dg ? DIAG : (k == hg ? HORZ : VERT);
        }
    }

    // Traceback to reconstruct optimal alignment
    j = 0;
    k = dn.at(m) <= vn.at(m) ? (dn.at(m) <= hn.at(m) ? DIAG : HORZ) : (vn.at(m) <= hn.at(m) ? VERT : HORZ);
    while(n>0 || m>0)
    {
        if(n < 0 || m < 0)
        {
        #ifdef DEBUG_MSA
            std::cerr << "something wrong here in newal(): n = " << n << ", m = " << m << std::endl;
        #endif
            break;
        }
        if(k == DIAG)
        {
            for(i = 1; i < I; ++i)
                NAL.at(i, j) = OAL.at(i, m);
            NAL.at(I, j++) = low + n;
            k = msa.dd.at(n--, m--);
        }
        else if(k == VERT)
        {
            for(i = 1; i < I; ++i)
                NAL.at(i, j) = 0;
            NAL.at(I, j++) = low + n;
            k = msa.vv.at(n--, m);
        }
        else
        {
            for(i = 1; i < I; ++i)
                NAL.at(i, j) = OAL.at(i, m);
            NAL.at(I, j++) = 0;
            k = msa.hh.at(n, m--);
        }
    }
    return j;
}

// Subroutine to calculate imposed cost for any pair of sequences
long long MSA_Ecalc::pcost(long long I, long long J, long long b)
{
    long long i;
    long long s = 0;
    for(i = 1; i <= b; ++i)
        s+= msa.D[ AL.at(I, i) ][ AL.at(J, i) ]+msa.T[ AL.at(I, i-1)!=DASH ][ AL.at(J, i-1)!=DASH][ AL.at(I, i)!=DASH][ AL.at(J, i)!=DASH ];
    i = 1;
    while( AL.at(I, i) == DASH && AL.at(J, i) == DASH)
        ++i;
    if(AL.at(I, i) == DASH || AL.at(J, i) == DASH)
        s+= (msa.GG-msa.G);
    i = b;
    while( AL.at(I, i) == DASH && AL.at(J, i) == DASH)
        --i;
    if(AL.at(I, i) == DASH || AL.at(J, i) == DASH)
        s += (msa.GG - msa.G);
    return s;
}


}// end of namespace msa_lib
