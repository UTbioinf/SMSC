#include "msa_bias.h"

namespace msa_lib
{
MSA_Bias::MSA_Bias(MSA& m): msa(m),
        Dij(2), B(2)
{}

void MSA_Bias::bias()
{
    vlist.reserve( (msa.S.size() - 1) << 1 );

    NODE* ancestor;
    double len;

    init_data(); // read in dist. matrix, initialize starting data structure
    minimize_Sij();  // determine opt node pairs, and make them neighbors
    // finish up turning initial data structure into binary tree by forming root node
    ancestor = makenode( ROOT, 0.0, node_index.at(0), node_index.at(1) );
    count2 = 1;
    len = dist(0, 1);
    len -= subdist(node_index.at(0), 0.0) / count2;
    count2 = 1;
    len -= subdist(node_index.at(1), 0.0) / count2;
    node_index.at(0)->ltop = len;

#ifdef DEBUG_MSA
    // print out rooted tree
    rpt(ancestor);
    std::cout << std::endl;
#endif

    // calculate weights
    whts();

    freevlist();
}

NODE* MSA_Bias::makenode(long long sqn, double ltop, NODE* lt, NODE* rt)
{
    NODE* ptr;

    ptr = new NODE();
    //vlist.at( vcount++ ) = ptr;
    vlist.push_back( ptr );     ++vcount;
    ptr->ltop = ltop;
    ptr->lt = lt;
    ptr->rt = rt;
    ptr->sqn = sqn;
    if(sqn < 0)
    {
        lt->bro = rt;
        rt->bro = lt;
        lt->par = rt->par = ptr;
    }
    return ptr;
}

double MSA_Bias::dist(long long i, long long j)
{
    count = 1;
    return ((double)rdist(node_index.at(i), node_index.at(j)) / count);
}

long long MSA_Bias::rdist(NODE* A, NODE* B)
{
    if(A->sqn < 0)
    {
        ++count;
        return (rdist(A->lt, B) + rdist(A->rt, B));
    }
    else if(B->sqn < 0)
    {
        ++count;
        return (rdist(A, B->lt) + rdist(A, B->rt));
    }
    return Dij.at(A->sqn, B->sqn);
}

double MSA_Bias::subdist(NODE* A, double total)
{
    if(A->sqn >= 0) return (total + A->ltop);
    ++count2;
    return (subdist(A->lt, A->ltop + total) + subdist(A->rt, A->ltop + total));
}

double MSA_Bias::br_len(long long i, long long j)
{
    double diz = 0, djz = 0;
    long long t;
    count2 = 1;
    for(t = 0; t < indexlen; ++t)
        if(t != i && t!= j)
        {
            diz += dist(i, t);
            djz += dist(j, t);
        }
    diz /= indexlen - 2;
    djz /= indexlen - 2;
    return ((dist(i, j) + diz - djz) / 2 - subdist(node_index.at(i), 0.0)/count2);
}

double MSA_Bias::compute_S(long long i, long long j)
{
    long long t, tt;
    double s1 = 0, s2 = 0;
    for(t = 0; t < indexlen; ++t)
        if(t != i && t != j)
            s1 += dist(i, t) + dist(j, t);
    s1 /= (indexlen - 2) << 1;
    for(t = 0; t < indexlen; ++t)
        for(tt = t + 1; tt < indexlen; ++tt)
            if(t != i && t != j && tt != i && tt != j)
                s2 += dist(t, tt);
    s2 /= indexlen -2;
    return (s1 + s2 + dist(i, j) / 2);
}

void MSA_Bias::coalesce(long long i, long long j)
{
    NODE* par;
    node_index.at(i)->ltop = br_len(i, j);
    node_index.at(j)->ltop = br_len(j, i);
    par = makenode(INOD, 0.0, node_index.at(i), node_index.at(j));
    node_index.at(i) = par;
    node_index.at(j) = node_index.at( indexlen - 1 );
    --indexlen;
}

void MSA_Bias::minimize_Sij()
{

    long long i, j, min_i = 0, min_j = 0;
    double tmp, min = BIG + 1;

    for(i = 0; i < indexlen; ++i)
        for(j = i + 1; j < indexlen; ++j)
        {
            tmp = compute_S(i, j);
            if(tmp < min)
            {
                min_i = i;
                min_j = j;
                min = tmp;
            }
        }
    coalesce(min_i, min_j);
    if(indexlen > 2)    minimize_Sij();
}

void MSA_Bias::init_data()
{
    long long i, j;
    indexlen = msa.S.size() - 1;
    Dij.set_dimensions( indexlen, indexlen);
    node_index.assign(indexlen, NULL);
    for(i = 0; i < indexlen -1; ++i)
        for(j = i + 1; j < indexlen; ++j)
            Dij.at( i, j ) = Dij.at(j, i) = msa.scale.at( j+1, i+1);
    for(i = vcount = 0; i < indexlen; ++i)
        node_index.at( i ) = makenode( i, 0.0, NULL, NULL );
}

#ifdef DEBUG_MSA
void MSA_Bias::rpt(NODE* A)
{
    if(A->sqn >= 0) std::cout << "Leaf #" << (1 + A->sqn) << "        Distance to parent = " << (A->ltop) << std::endl;
    else
    {
        if(A->sqn == ROOT)
            std::cout << "---------------- Tree given from ancestor ----------------" << std::endl;
        else
            std::cout << "Internal Node  Distance to parent = " << (A->ltop) << std::endl;
        std::cout << "On the left:   ";
        rpt(A->lt);
        std::cout << "On the right:  ";
        rpt(A->rt);
    }
}
#endif

void MSA_Bias::whts()
{
    long long i, j;
    NODE* no;
    double sm;

    size_t K = msa.S.size() - 1;
    B.set_dimensions(K, K);
    for(size_t pN = 0; (no = vlist.at( pN ))->sqn > INOD; ++pN)
        trace( vlist.at(pN), 1.0, no->par, no->bro );
    for(sm = BIG, j = 1; j < K; ++j)
        for(i = 0; i < j; ++i)
            if(B.at( i, j ) < sm)
                sm = B.at(i, j);
    for(i = 0; i < K -1; ++i)
        for(j = i + 1; j < K; ++j)
            msa.scale.at(i+1, j + 1) = B.at(i, j) / sm + 0.5;
}

void MSA_Bias::trace(NODE* pN, double prod, NODE* no, NODE* sis)
{
    if(no->sqn > INOD)  B.at(pN->sqn, no->sqn) = prod;
    else if(sis == NULL)
    {
        trace(pN, prod/2, no->rt, NULL);
        trace(pN, prod/2, no->lt, NULL);
    }
    else if(no->sqn != ROOT)
    {
        trace(pN, prod/2, sis, NULL);
        trace(pN, prod/2, no->par, no->bro);
    }
    else    trace(pN, prod, sis, NULL);
}

void MSA_Bias::freevlist()
{
    while(vcount>0)
        delete vlist.at( --vcount );
    vlist.clear();
}


}// end namespace msa_lib
