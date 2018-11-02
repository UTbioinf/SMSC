#ifndef _MSA_ECALC_H
#define _MSA_ECALC_H

#include "msa_base.h"
#include "msa_lib.h"

namespace msa_lib
{
class MSA_Ecalc
{
private:
    MSA& msa;
    int K;

    std::vector<const std::string*> SS;    // [NUMBER+1]
    NDArray<long long> AL;                // [NUMBER+1][2*LENGTH]
    NDArray<long long> NAL;               // [NUMBER+1][2*LENGTH]
    NDArray<long long> OAL;               // [NUMBER+1][2*LENGTH]
    std::vector<long long> Ord;           // [NUMBER+1]
#ifdef DEBUG_MSA
    std::string A;
#else
    const char* A;
#endif
public:
    MSA_Ecalc(MSA& m);
    void ecalc(long long len, bool optimal);
    void order(long long oldp1, long long p1);
    long long align(long long oldp1, long long p1, long long len);
    long long newal(long long I, long long low, long long n, long long m);
    long long pcost(long long I, long long J, long long b);    
};

} // end of namespace msa_lib
#endif
