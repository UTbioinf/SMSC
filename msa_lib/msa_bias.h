#ifndef _MSA_BIAS_H
#define _MSA_BIAS_H

#include "msa_base.h"
#include "msa_lib.h"

namespace msa_lib
{
class MSA_Bias
{
private:
    MSA& msa;

    std::vector<NODE* > node_index;
    NDArray<long long> Dij;
    long long vcount, indexlen, count, count2;
    std::vector<NODE *> vlist;
    NDArray<double> B;
public:
    MSA_Bias(MSA& m);
    void bias();
    NODE* makenode(long long sqn, double ltop, NODE* lt, NODE* rt);
    double dist(long long i, long long j);
    long long rdist(NODE* A, NODE* B);
    double subdist(NODE* A, double total);
    double br_len(long long i, long long j);
    double compute_S(long long i, long long j);
    void coalesce(long long i, long long j);
    void minimize_Sij();
    void init_data();
#ifdef DEBUG_MSA
    void rpt(NODE *A);
#endif
    void whts();
    void trace(NODE* pN, double prod, NODE* no, NODE* sis);
    void freevlist();    
};
}// end namespace msa_lib

#endif
