#ifndef _DELTA_FILTER_H
#define _DELTA_FILTER_H

#include <string>

#include "delta.h"
#include "tigrinc.h"
//#include "postnuc.h"

#define DELTA_FILTER_SET_e 1
#define DELTA_FILTER_SET_g 2
#define DELTA_FILTER_SET_i 4
#define DELTA_FILTER_SET_l 8
#define DELTA_FILTER_SET_o 0x10
#define DELTA_FILTER_SET_q 0x20
#define DELTA_FILTER_SET_r 0x40
#define DELTA_FILTER_SET_u 0x80
#define DELTA_FILTER_SET_m 0x100
#define DELTA_FILTER_SET_1 0x200

typedef struct
{
    float e;            // negligible alignment score
    float i;            // minimum identity
    float l;            // minimum alignment length
    float o;            // maximum overlap as % of align len
    float u;            // minimum unique
    
    int options;
}delta_filter_Params;


typedef struct
{
    void* res;  // store result in this place
    void (*store_header)(void*, const char*, const char*, long int, long int);
    void (*new_cluster)(void*, long int, long int, long int, long int, long int, long int, long int);
    void (*add_delta)(void*, long int);
}delta_filter_ProcessResult;

class DeltaFilter: public tigrinc::DeltaGraph_t
{
private:
    std::string         OPT_AlignName;// alignment name parameter
    
    bool           OPT_QLIS;     // do query based LIS
    bool           OPT_RLIS;     // do reference based LIS
    bool           OPT_GLIS;     // do global LIS
    bool           OPT_1to1;     // do 1-to-1 alignment
    bool           OPT_MtoM;     // do M-to-M alignment
    long int       OPT_MinLength;    // minimum alignment length
    float          OPT_MinIdentity;  // minimum %identity
    float          OPT_MinUnique;    // minimum %unique
    float          OPT_MaxOverlap;   // maximum olap as % of align len
    float          OPT_Epsilon;      // negligible alignment score

    bool           is_first_record;

    tigrinc::DeltaRecord_t record_m;
    float record_m_total;

    void insert_dep();
public:
    DeltaFilter();
    void set_parameters(const delta_filter_Params& params);
    void reset_parameters();
    void set_refpath(const std::string& refpath);
    void set_qrypath(const std::string& qrypath);
    void set_datatype(const std::string& datatype);
    void add_record_header(const std::string& idR, const std::string& idQ, long lenR, long lenQ);
    void add_alignment(long sR, long eR, long sQ, long eQ, long idyc, long simc, long stpc, const std::vector<long>& all_deltas);
    void add_alignment(long sR, long eR, long sQ, long eQ, long idyc, long simc, long stpc);
    void add_indel(long d);
    void add_finished(); // call this only if you called the second add_alignment() and the add_indel();
    void run_delta_filter();
    void get_result(delta_filter_ProcessResult& processFilter);
    void clear();
};

#endif
