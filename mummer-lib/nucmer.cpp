#include <cassert>

#include "nucmer.h"

namespace tigrinc
{
static Suffixtree stree;
static Match_t mt;
static nucmer_Params nuc_params;
static Uint ref_offset = 0;
static Uint ref_cut_L = -1;
static Uint ref_cut_R = -1;

/* =========== default callback functions ==============*/
static Sint save_maximalmatch_rightorder(void *info, Uint matchlength, Uint refstart, Uint seqnum, Uint qrystart)
{
    if(ref_cut_L == -1)
    {
        vector<Match_t> *triples = (vector<Match_t> *)(info);
        mt.Start1 = refstart + 1 + ref_offset;
        mt.Start2 = qrystart + 1;
        mt.Len = matchlength;

        triples->push_back( mt );
    }
    else
    {
        if(refstart > ref_cut_R || refstart + matchlength <= ref_cut_L)
        {
            vector<Match_t> *triples = (vector<Match_t> *)(info);
            mt.Start1 = refstart + 1 + ref_offset;
            mt.Start2 = qrystart + 1;
            mt.Len = matchlength;

            triples->push_back( mt );
        }
    }
    return 0;
}

static Sint save_maximalmatch_reverseorder(void *info, Uint matchlength, Uint qrystart, Uint seqnum, Uint refstart)
{
    if(ref_cut_L == -1)
    {
        vector<Match_t> *triples = (vector<Match_t> *)(info);
        mt.Start1 = refstart + 1 + ref_offset;
        mt.Start2 = qrystart + 1;
        mt.Len = matchlength;

        triples->push_back( mt );
    }
    else
    {
        if(refstart > ref_cut_R || refstart + matchlength <= ref_cut_L)
        {
            vector<Match_t> *triples = (vector<Match_t> *)(info);
            mt.Start1 = refstart + 1 + ref_offset;
            mt.Start2 = qrystart + 1;
            mt.Len = matchlength;

            triples->push_back( mt );
        }
    }
    return 0;
}

static void add_mgaps_item(void* c, int new_group, long int Start1, long int Start2, long int Len,
        long int adj, long int s1_olap, long int s2_olap)
{
    if(new_group)
        postnuc_add_Match(Start1, Start2, Len, true);
    else
        postnuc_add_Match(Start1, Start2, Len);
}


static void save_header_postnuc(void* res, const char* ref_header, const char* qry_header, 
        long int ref_len, long int qry_len)
{
    DeltaFilter* df = (DeltaFilter*) res;
    df->add_record_header( ref_header, qry_header, ref_len, qry_len );
}

static void new_cluster_postnuc(void* res, long int sA, long int eA, 
        long int sB, long int eB, 
        long int u1, long int u2, long int u3)
{
    DeltaFilter* df = (DeltaFilter *) res;
    df->add_alignment( sA, eA, sB, eB, u1, u2, u3 );
}

static void add_delta_postnuc(void* res, long int d)
{
    DeltaFilter* df = (DeltaFilter *) res;
    df->add_indel( d );
}

static void save_header_deltafilter(void* res, const char* ref_header, const char* qry_header, 
        long int ref_len, long int qry_len)
{
    Nucmer_Deltas* d = (Nucmer_Deltas* )res;
    d->push_back( Nucmer_Delta() );
    Nucmer_Delta& t_d = d->back();

    t_d.ref_tag = ref_header;
    t_d.qry_tag = qry_header;
    t_d.ref_len = ref_len;
    t_d.qry_len = qry_len;
}

static void new_cluster_deltafilter(void* res, long int sA, long int eA, 
        long int sB, long int eB, 
        long int u1, long int u2, long int u3)
{
    Nucmer_Deltas* d = (Nucmer_Deltas* )res;

    d->back().clusters.push_back( Nucmer_Cluster() );
    Nucmer_Cluster& clu = d->back().clusters.back();

    clu.ref_start = sA;
    clu.ref_end = eA;
    clu.qry_start = sB;
    clu.qry_end = eB;
    clu.unknown[0] = u1;
    clu.unknown[1] = u2;
    clu.unknown[2] = u3;
}

static void add_delta_deltafilter(void* res, long int d)
{
    Nucmer_Deltas* dd = (Nucmer_Deltas*) res;
    dd->back().clusters.back().indels.push_back( d );
}

/* ============ classes ===============*/

Nucmer_Cluster::Nucmer_Cluster():
    ref_start(0), ref_end(0),
    qry_start(0), qry_end(0),
    is_alive(true)
{
}

void Nucmer_Cluster::swap(Nucmer_Cluster& clu)
{
    
    if(&clu == this)
        return;
    ::swap(this->is_alive, clu.is_alive);
    ::swap(this->ref_start, clu.ref_start);
    ::swap(this->ref_end, clu.ref_end);
    ::swap(this->qry_start, clu.qry_start);
    ::swap(this->qry_end, clu.qry_end);
    ::swap(this->unknown[0], clu.unknown[0]);
    ::swap(this->unknown[1], clu.unknown[1]);
    ::swap(this->unknown[2], clu.unknown[2]);
    this->indels.swap( clu.indels );
}

void Nucmer_Deltas::print_deltas(ostream& out)
{
    for(vector<Nucmer_Delta>::iterator d_it = this->begin(); d_it != this->end(); ++d_it)
    {
        out << '>' << (d_it->ref_tag) << ' ' << (d_it->qry_tag) << ' ' << (d_it->ref_len) << ' ' << (d_it->qry_len) << endl;
        for(vector<Nucmer_Cluster>::iterator c_it = d_it->clusters.begin(); c_it != d_it->clusters.end(); ++c_it)
        {
            out << (c_it->ref_start) << ' '
                << (c_it->ref_end) << ' '
                << (c_it->qry_start) << ' '
                << (c_it->qry_end) << endl;
            for(vector<long>::iterator id_it = c_it->indels.begin(); id_it != c_it->indels.end(); ++id_it)
                out << (*id_it) << endl;
            out << 0 << endl;
        }
    }
}

bool Nucmer_Deltas::trim_left(Nucmer_Cluster& cluster_left, Nucmer_Cluster& cluster_right)
{
    if(!cluster_left.is_alive)
        return true; // continue trim

    if(cluster_left.ref_end < cluster_right.ref_start && cluster_left.qry_end < cluster_right.qry_start)
        return false;

    if(cluster_left.ref_start >= cluster_right.ref_start || cluster_left.qry_start >= cluster_right.qry_start)
    {
        cluster_left.is_alive = false;
        cluster_left.indels.clear();
        return true;
    }

    if(cluster_right.ref_end <= cluster_left.ref_end || cluster_right.qry_end <= cluster_left.qry_end)
    {
        cluster_right.is_alive =false;
        cluster_right.indels.clear();
        return false;
    }

    // remove redundant indels
    long t_ref = cluster_left.ref_start - 1;
    long t_qry = cluster_left.qry_start - 1;
    for(vector<long>::iterator d_it = cluster_left.indels.begin();
            d_it != cluster_left.indels.end(); ++d_it)
    {
        if(*d_it > 0)
        {
            t_ref += *d_it;
            t_qry += *d_it - 1;
        }
        else
        {
            t_ref += -(*d_it) - 1;
            t_qry += -(*d_it);
        }

        if(t_ref < cluster_right.ref_start && t_qry < cluster_right.qry_start)
            continue;

        if(*d_it > 0)
            --t_ref;
        else
            --t_qry;

        if(t_ref >= cluster_right.ref_start || t_qry >= cluster_right.qry_start)
        {
            long max_trim = ::max(t_ref - cluster_right.ref_start + 1, t_qry - cluster_right.qry_start + 1);
            cluster_left.ref_end = t_ref - max_trim;
            cluster_left.qry_end = t_qry - max_trim;
        }
        else
        {
            cluster_left.ref_end = t_ref;
            cluster_left.qry_end = t_qry;
        }

        if(cluster_left.ref_end <= cluster_left.ref_start || cluster_left.qry_end <= cluster_left.qry_start)
        {
            cluster_left.is_alive = false;
            cluster_left.indels.clear();
        }
        else
        {
            cluster_left.indels.erase(d_it, cluster_left.indels.end());
        }

        return true;
    }

    if(cluster_left.ref_end >= cluster_right.ref_start || cluster_left.qry_end >= cluster_right.qry_start)
    {
        long max_trim = ::max(cluster_left.ref_end - cluster_right.ref_start + 1,
                                cluster_left.qry_end - cluster_right.qry_start + 1);
        cluster_left.ref_end -= max_trim;
        cluster_left.qry_end -= max_trim;

        if(cluster_left.ref_end <= cluster_left.ref_start || cluster_left.qry_end <= cluster_left.qry_start)
        {
            cluster_left.is_alive = false;
            cluster_left.indels.clear();
        }
    }
    return true;
}

void Nucmer_Deltas::trim()
{
    for(vector<Nucmer_Delta>::iterator delta_it = this->begin(); delta_it != this->end(); ++delta_it)
    {
        vector<Nucmer_Cluster>& t_clusters = delta_it -> clusters;
        if(t_clusters.size() < 2)
            continue;

        for(size_t i = 1; i < t_clusters.size(); ++i)
        {
            size_t j = i - 1;
            bool continue_trim = false;
            do
            {
                continue_trim = trim_left(t_clusters[j], t_clusters[i]);
                if(j == 0)
                    break;
                else
                    --j;
            }while(continue_trim);
        }
    }
    remove_dead();
}

void Nucmer_Deltas::remove_dead()
{
    for(vector<Nucmer_Delta>::iterator d_it = this->begin(); d_it != this->end(); ++d_it)
    {
        vector<Nucmer_Cluster>::iterator cur_it = d_it->clusters.begin();
        vector<Nucmer_Cluster>::iterator c_it = d_it->clusters.begin();

        for(; c_it != d_it->clusters.end(); ++c_it)
        {
            if(c_it->is_alive)
            {
                c_it->swap( *cur_it );
                ++cur_it;
            }
        }

        d_it->clusters.erase(cur_it, d_it->clusters.end());
    }
}

void Nucmer_Deltas::check_deltas()
{
    for(vector<Nucmer_Delta>::const_iterator d_it = this->begin(); d_it != this->end(); ++d_it)
    {
        for(vector<Nucmer_Cluster>::const_iterator clu_it = d_it->clusters.begin();
                clu_it != d_it->clusters.end(); ++clu_it)
        {
            if(!clu_it->is_alive)
                continue;
            long ref_pos = -1;
            long qry_pos = -1;
            for(vector<long>::const_iterator indel_it = clu_it->indels.begin();
                    indel_it != clu_it->indels.end(); ++indel_it)
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

            assert(clu_it->ref_end - clu_it->ref_start - ref_pos == clu_it->qry_end - clu_it->qry_start - qry_pos);
            assert(clu_it->ref_end - clu_it->ref_start - ref_pos >= 0);
        }
    }
}

void compute_maxmat(Uchar* ref, Uchar* qry, Uint ref_len, Uint qry_len, vector<Match_t>* triples)
{
    if(nuc_params.save_maximalmatch)
    {
        if(constructstree(&stree, ref, ref_len) != 0)
        {
            cerr << "Construct suffix tree fails" << endl;
            exit(1);
        }

        if(findmaxmatches(&stree,
                nuc_params.minmatchlength,
                nuc_params.save_maximalmatch,
                (void*)triples,
                qry,
                qry_len,
                0) !=0 )
        {
            cerr << "Find maxmatches error" << endl;
            exit(1);
        }
        freestree(&stree);
    }
    else
    {
        if(ref_len < qry_len)
        {
            if(constructstree(&stree, ref, ref_len) != 0)
            {
                cerr << "Construct suffix tree fails" << endl;
                exit(1);
            }

            if(findmaxmatches(&stree,
                    nuc_params.minmatchlength,
                    save_maximalmatch_rightorder,
                    (void *)triples,
                    qry,
                    qry_len,
                    0) != 0)
            {
                cerr << "Find maxmatches error" << endl;
                exit(1);
            }
        }
        else
        {
            if(constructstree(&stree, qry, qry_len) != 0)
            {
                cerr << "Construct suffix tree fails" << endl;
                exit(1);
            }
            if(findmaxmatches(&stree,
                    nuc_params.minmatchlength,
                    save_maximalmatch_reverseorder,
                    (void *)triples,
                    ref,
                    ref_len,
                    0) != 0)
            {
                cerr << "Find maxmatches error" << endl;
                exit(1);
            }
        }
        freestree(&stree);
    }
}

void compute_mgaps(vector<Match_t>& triples)
{
    mgaps_Params params={0};
    if(nuc_params.mgp_params)
        params = *nuc_params.mgp_params;
    else
    {
        params.l = 50;
        params.s = 20;
        params.d = 5;
        params.f = .12;
        params.options = MGAPS_SET_L | MGAPS_SET_S | MGAPS_SET_D | MGAPS_SET_F;
    }
    mgaps_set_parameters(&params);

    if(nuc_params.add_mgaps_item)
        mgaps_Process_Matches(&triples.front(), triples.size()-1, "lol", nuc_params.add_mgaps_item, NULL);
    else
        mgaps_Process_Matches(&triples.front(), triples.size()-1, "lol", add_mgaps_item, NULL);
}

void compute_deltas(Uchar* ref, Uchar* qry, Uint ref_len, Uint qry_len, Nucmer_Deltas& deltas, bool do_trim)
{
    vector<Match_t> triples;
    triples.push_back(mt);

    DeltaFilter df;
    delta_filter_Params df_params;
    delta_filter_ProcessResult processFilter;
    if(nuc_params.dfp_params)
        df_params = *nuc_params.dfp_params;
    else
        df_params.options = DELTA_FILTER_SET_g;
    if(nuc_params.dfp_processFilter)
        processFilter = *nuc_params.dfp_processFilter;
    else
    {
        processFilter.res = &deltas;
        processFilter.store_header = save_header_deltafilter;
        processFilter.new_cluster = new_cluster_deltafilter;
        processFilter.add_delta = add_delta_deltafilter;
    }
    df.set_datatype("NUCMER");
    df.set_parameters(df_params);

    postnuc_ProcessResult deltaProcess;
    postnuc_Params params;
    if(nuc_params.dtp_process)
        deltaProcess = *nuc_params.dtp_process;
    else
    {
        deltaProcess.res = &df;
        deltaProcess.store_header = save_header_postnuc;
        deltaProcess.new_cluster = new_cluster_postnuc;
        deltaProcess.add_delta = add_delta_postnuc;
    }
    if(nuc_params.pnp_params)
        params = *nuc_params.pnp_params;
    else
    {
        params.b = 150;
        params.options = POSTNUC_SET_b;
    }
    postnuc_init(&params);
    postnuc_set_ref("ref", (long int)ref_len, (char *)ref);
    postnuc_set_qry("qry", (long int)qry_len, (char *)qry);

    compute_maxmat(ref, qry, ref_len, qry_len, &triples);
    compute_mgaps(triples);
    postnuc_computeDelta( deltaProcess, true );
    df.add_finished();

    df.run_delta_filter();

    df.get_result(processFilter);

    if(do_trim)
    {
        deltas.trim();
    }
    
}

void compute_deltas_once(Suffixtree& stree, Uchar* ref, Uchar* qry, Uint ref_len, Uint qry_len, Nucmer_Deltas& deltas, bool do_trim)
{
    vector<Match_t> triples;
    triples.push_back(mt);

    DeltaFilter df;
    delta_filter_Params df_params;
    delta_filter_ProcessResult processFilter;
    if(nuc_params.dfp_params)
        df_params = *nuc_params.dfp_params;
    else
        df_params.options = DELTA_FILTER_SET_g;
    if(nuc_params.dfp_processFilter)
        processFilter = *nuc_params.dfp_processFilter;
    else
    {
        processFilter.res = &deltas;
        processFilter.store_header = save_header_deltafilter;
        processFilter.new_cluster = new_cluster_deltafilter;
        processFilter.add_delta = add_delta_deltafilter;
    }
    df.set_datatype("NUCMER");
    df.set_parameters(df_params);

    postnuc_ProcessResult deltaProcess;
    postnuc_Params params;
    if(nuc_params.dtp_process)
        deltaProcess = *nuc_params.dtp_process;
    else
    {
        deltaProcess.res = &df;
        deltaProcess.store_header = save_header_postnuc;
        deltaProcess.new_cluster = new_cluster_postnuc;
        deltaProcess.add_delta = add_delta_postnuc;
    }
    if(nuc_params.pnp_params)
        params = *nuc_params.pnp_params;
    else
    {
        params.b = 150;
        params.options = POSTNUC_SET_b;
    }
    postnuc_init(&params);
    postnuc_set_ref("ref", (long int)ref_len, (char *)ref);
    postnuc_set_qry("qry", (long int)qry_len, (char *)qry);

    // compute maxmat
    if(nuc_params.save_maximalmatch)
    {
        if(findmaxmatches(&stree,
                nuc_params.minmatchlength,
                nuc_params.save_maximalmatch,
                (void*)&triples,
                qry,
                qry_len,
                0) !=0 )
        {
            cerr << "Find maxmatches error" << endl;
            exit(1);
        }
    }
    else
    {
        if(findmaxmatches(&stree,
                    nuc_params.minmatchlength,
                    save_maximalmatch_rightorder,
                    (void *)&triples,
                    qry,
                    qry_len,
                    0) != 0)
        {
            cerr << "Find maxmatches error" << endl;
            exit(1);
        }
    }

    compute_mgaps(triples);
    postnuc_computeDelta( deltaProcess, true );
    df.add_finished();

    df.run_delta_filter();

    df.get_result(processFilter);

    if(do_trim)
    {
        deltas.trim();
    }
}

}

/* =============== out of namespace tigrinc =================*/

nucmer_Params::nucmer_Params()
{
    init();
}

void nucmer_Params::init()
{
    dfp_params = NULL;
    dfp_processFilter = NULL;
    dtp_process = NULL;
    pnp_params = NULL;
    minmatchlength = 20;
    save_maximalmatch = NULL;
    mgp_params = NULL;
    add_mgaps_item = NULL;
    
}

void nucmer_set_parameters(nucmer_Params& p)
{
    tigrinc::nuc_params = p;
}

void nucmer_reset_parameters()
{
    tigrinc::nuc_params.init();
}

void nucmer_compute_deltas(const string& ref_str, const string& qry_str, nucmer_Deltas& deltas, bool do_trim/* = false*/)
{
    Uchar* ref = reinterpret_cast<Uchar *>(const_cast<char *>(ref_str.c_str()));
    Uchar* qry = reinterpret_cast<Uchar *>(const_cast<char *>(qry_str.c_str()));
    Uint ref_len = ref_str.length();
    Uint qry_len = qry_str.length();

    tigrinc::compute_deltas(ref, qry, ref_len, qry_len, deltas, do_trim);
}

void nucmer_compute_deltas_once(Suffixtree& stree, Uchar* ref, Uint ref_len, const string& qry_str, nucmer_Deltas& deltas, bool do_trim/* = false*/)
{
    Uchar* qry = reinterpret_cast<Uchar *>(const_cast<char *>(qry_str.c_str()));
    Uint qry_len = qry_str.length();

    tigrinc::compute_deltas_once(stree, ref, qry, ref_len, qry_len, deltas, do_trim);
}

void nucmer_alignall(const string& ref_str, const string& qry_str, nucmer_Deltas& deltas)
{
    nucmer_compute_deltas(ref_str, qry_str, deltas, true);
}

void nucmer_alignleft(const string& ref_str, const string& qry_str, nucmer_Deltas& deltas)
{
    Uint ref_len = min(ref_str.length(), qry_str.length());
    Uint qry_len = qry_str.length();

    Uchar* ref = new Uchar[ ref_len + 1];
    Uchar* qry = reinterpret_cast<Uchar *>(const_cast<char *>(qry_str.c_str()));

    for(Uint i = 0; i < ref_len; ++i)
        ref[i] = ref_str[i];
    ref[ref_len] = 0;

    tigrinc::compute_deltas(ref, qry, ref_len, qry_len, deltas, true);

    delete[] ref;
}

void nucmer_alignright(const string& ref_str, const string& qry_str, nucmer_Deltas& deltas)
{
    
    Uint ref_len = min(ref_str.length(), qry_str.length());
    Uint qry_len = qry_str.length();
    if(ref_len < ref_str.length())
        tigrinc::ref_offset = ref_str.length() - ref_len;

    Uchar* ref = new Uchar[ ref_len + 1];
    Uchar* qry = reinterpret_cast<Uchar *>(const_cast<char *>(qry_str.c_str()));

    for(Uint i = 0; i < ref_len; ++i)
        ref[i] = ref_str[i + tigrinc::ref_offset];
    ref[ref_len] = 0;

    tigrinc::compute_deltas(ref, qry, ref_len, qry_len, deltas, true);

    delete[] ref;
    tigrinc::ref_offset = 0;
}

void nucmer_alignends(const string& ref_str, const string& qry_str, nucmer_Deltas& deltas)
{
    Uchar* ref = reinterpret_cast<Uchar *>(const_cast<char *>(ref_str.c_str()));
    Uchar* qry = reinterpret_cast<Uchar *>(const_cast<char *>(qry_str.c_str()));
    Uint ref_len = ref_str.length();
    Uint qry_len = qry_str.length();

    if(ref_len < qry_len + qry_len)
        tigrinc::compute_deltas(ref, qry, ref_len, qry_len, deltas, true);
    else
    {
        tigrinc::ref_cut_L = qry_len;
        tigrinc::ref_cut_R = ref_len - qry_len - 1;

        tigrinc::compute_deltas(ref, qry, ref_len, qry_len, deltas, true);

        tigrinc::ref_cut_L = -1;
        tigrinc::ref_cut_R = -1;
    }
}

void nucmer_alignends_once(Suffixtree& stree, Uchar* ref, Uint ref_len, const string& qry_str, nucmer_Deltas& deltas)
{
    Uchar* qry = reinterpret_cast<Uchar *>(const_cast<char *>(qry_str.c_str()));
    Uint qry_len = qry_str.length();

    if(ref_len < qry_len + qry_len)
        tigrinc::compute_deltas_once(stree, ref, qry, ref_len, qry_len, deltas, true);
    else
    {
        tigrinc::ref_cut_L = qry_len;
        tigrinc::ref_cut_R = ref_len - qry_len - 1;

        tigrinc::compute_deltas_once(stree, ref, qry, ref_len, qry_len, deltas, true);

        tigrinc::ref_cut_L = -1;
        tigrinc::ref_cut_R = -1;
    }
}
