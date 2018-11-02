#ifndef _NUCMER_H
#define _NUCMER_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include "types.h"
#include "debugdef.h"
#include "errordef.h"
#include "protodef.h"
#include "spacedef.h"
#include "streedef.h"
#include "maxmatdef.h"
#include "mgaps.h"
#include "postnuc.h"
#include "delta-filter.h"

using namespace std;

namespace tigrinc
{
typedef Sint (*Save_Maximalmatch)(void*, Uint, Uint, Uint, Uint);

class Nucmer_Cluster
{
public:
    bool is_alive;
    long long ref_start, ref_end;
    long long qry_start, qry_end;
    long long unknown[3];
    vector<long> indels;
public:
    Nucmer_Cluster();
    void swap(Nucmer_Cluster& clu);
};

class Nucmer_Delta
{
public:
    string ref_tag, qry_tag;
    long long ref_len, qry_len;
    vector<Nucmer_Cluster> clusters;
};

class Nucmer_Deltas: public vector<Nucmer_Delta>
{
public:
    void print_deltas(ostream& out);
    bool trim_left(Nucmer_Cluster& cluster_left, Nucmer_Cluster& cluster_right);
    void trim();
    void remove_dead();
    void check_deltas();
};

void compute_maxmat(Uchar* ref, Uchar* qry, Uint ref_len, Uint qry_len, vector<Match_t>* triples);
void compute_mgaps(vector<Match_t>& triples);
void compute_deltas(Uchar* ref, Uchar* qry, Uint ref_len, Uint qry_len, Nucmer_Deltas& deltas, bool do_trim);
void compute_deltas_once(Suffixtree& stree, Uchar* ref, Uchar* qry, Uint ref_len, Uint qry_len, Nucmer_Deltas& deltas, bool do_trim);
}

/* out of namespace tigrinc */
typedef tigrinc::Nucmer_Deltas nucmer_Deltas;
class nucmer_Params
{
public:
    delta_filter_Params* dfp_params;
    delta_filter_ProcessResult* dfp_processFilter;
    postnuc_ProcessResult* dtp_process;
    postnuc_Params* pnp_params;
    Uint minmatchlength;
    tigrinc::Save_Maximalmatch save_maximalmatch;
    mgaps_Params *mgp_params;
    void (*add_mgaps_item)(void*, int, long int, long int, long int,
            long int, long int, long int);
public:
    nucmer_Params();
    void init();
};

void nucmer_set_parameters(nucmer_Params& p);
void nucmer_reset_parameters();
void nucmer_compute_deltas(const string& ref_str, const string& qry_str, nucmer_Deltas& deltas, bool do_trim=false);
void nucmer_compute_deltas_once(Suffixtree& stree, Uchar* ref, Uint ref_len, const string& qry_str, nucmer_Deltas& deltas, bool do_trim=false);
void nucmer_alignall(const string& ref_str, const string& qry_str, nucmer_Deltas& deltas);
void nucmer_alignleft(const string& ref_str, const string& qry_str, nucmer_Deltas& deltas);
void nucmer_alignright(const string& ref_str, const string& qry_str, nucmer_Deltas& deltas);
void nucmer_alignends(const string& ref_str, const string& qry_str, nucmer_Deltas& deltas);
void nucmer_alignends_once(Suffixtree& stree, Uchar* ref, Uint ref_len, const string& qry_str, nucmer_Deltas& deltas);

#endif
