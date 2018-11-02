//------------------------------------------------------------------------------
//   Modified by: zijuexiansheng
//   Programmer: Adam M Phillippy, The Institute for Genomic Research
//         File: delta-filter.cc
//         Date: 09 / 22 / 2004
//
//        Usage: delta-filter  [options]  <deltafile>
//               Try 'show-coords -h' for more information
//
//  Description: For use in conjunction with the MUMmer package.
//              "delta-filter" cleans delta alignment files by filtering
//             alignments that fail to meet the specifications required by the
//            command line switches. Produces a new, filtered delta file on
//           stdout, and works for both nucleotide and amino-acid alignments.
//
//------------------------------------------------------------------------------

#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>

#include "delta-filter.h"

using namespace std;

void DeltaFilter::insert_dep()
{
    tigrinc::DeltaEdge_t* dep = new tigrinc::DeltaEdge_t();
    pair<map<string, tigrinc::DeltaNode_t>::iterator, bool> insret;

    //-- Find the reference node in the graph, add a new one if necessary
    insret = refnodes.insert(map<string, tigrinc::DeltaNode_t>::value_type
            (record_m.idR, tigrinc::DeltaNode_t()));
    dep->refnode = &((insret.first)->second);
    //-- If a new reference node
    if ( insret.second )
    {
        dep->refnode->id  = &((insret.first)->first);
        dep->refnode->len = record_m.lenR;
    }


    //-- Find the query node in the graph, add a new one if necessary
    insret = qrynodes.insert(map<string, tigrinc::DeltaNode_t>::value_type
            (record_m.idQ, tigrinc::DeltaNode_t()));
    dep->qrynode = &((insret.first)->second);
    //-- If a new query node
    if ( insret.second )
    {
        dep->qrynode->id  = &((insret.first)->first);
        dep->qrynode->len = record_m.lenQ;
    }


    //-- Build the edge
    dep->build(record_m);
    dep->refnode->edges.push_back(dep);
    dep->qrynode->edges.push_back(dep);
    
    record_m.clear();
}

// Definition of DeltaFilter
DeltaFilter::DeltaFilter():
    OPT_QLIS(false),     // do query based LIS
    OPT_RLIS(false),     // do reference based LIS
    OPT_GLIS(false),     // do global LIS
    OPT_1to1(false),     // do 1-to-1 alignment
    OPT_MtoM(false),     // do M-to-M alignment
    OPT_MinLength(0),         // minimum alignment length
    OPT_MinIdentity(0.0),     // minimum %identity
    OPT_MinUnique(0.0),       // minimum %unique
    OPT_MaxOverlap(100.0),    // maximum olap as % of align len
    OPT_Epsilon(-1.0),        // negligible alignment score
    is_first_record(true)
{
    srand(1);
}


void DeltaFilter::set_parameters(const delta_filter_Params& params)
{
    if(params.options & DELTA_FILTER_SET_e)
        OPT_Epsilon = params.e;
    if(params.options & DELTA_FILTER_SET_g)
        OPT_GLIS = true;
    if(params.options & DELTA_FILTER_SET_i)
        OPT_MinIdentity = params.i;
    if(params.options & DELTA_FILTER_SET_l)
        OPT_MinLength = params.l;
    if(params.options & DELTA_FILTER_SET_o)
        OPT_MaxOverlap = params.o;
    if(params.options & DELTA_FILTER_SET_q)
        OPT_QLIS = true;
    if(params.options & DELTA_FILTER_SET_r)
        OPT_RLIS = true;
    if(params.options & DELTA_FILTER_SET_u)
        OPT_MinUnique = params.u;
    if(params.options & DELTA_FILTER_SET_m)
        OPT_MtoM = true;
    if(params.options & DELTA_FILTER_SET_1)
        OPT_1to1 = true;

    if ( OPT_MinIdentity < 0.0  ||  OPT_MinIdentity > 100.0 )
    {
        cerr << "ERROR: Minimum identity must be within the range [0, 100]\n";
        exit(1);
    }

    if ( OPT_MinLength < 0 )
    {
        cerr << "ERROR: Minimum length must be greater than or equal to zero\n";
        exit(1);
    }

    if ( OPT_MinUnique < 0.0  ||  OPT_MinUnique > 100.0 )
    {
        cerr << "ERROR: Minimum uniqueness must be within the range [0, 100]\n";
        exit(1);
    }

    if ( OPT_MaxOverlap < 0.0  ||  OPT_MaxOverlap > 100.0 )
    {
        cerr << "ERROR: Maximum overlap must be within the range [0, 100]\n";
        exit(1);
    }
}

void DeltaFilter::reset_parameters()
{
    OPT_QLIS = false;     // do query based LIS
    OPT_RLIS = false;     // do reference based LIS
    OPT_GLIS = false;     // do global LIS
    OPT_1to1 = false;     // do 1-to-1 alignment
    OPT_MtoM = false;     // do M-to-M alignment
    OPT_MinLength = 0;         // minimum alignment length
    OPT_MinIdentity = 0.0;     // minimum %identity
    OPT_MinUnique = 0.0;       // minimum %unique
    OPT_MaxOverlap = 100.0;    // maximum olap as % of align len
    OPT_Epsilon = -1.0;      // negligible alignment score
}

void DeltaFilter::set_refpath(const string& refpath)
{
    this->refpath = refpath;
}

void DeltaFilter::set_qrypath(const string& qrypath)
{
    this->qrypath = qrypath;
}

void DeltaFilter::set_datatype(const string& datatype)
{
    if(datatype == tigrinc::NUCMER_STRING)
        this->datatype = tigrinc::NUCMER_DATA;
    else if(datatype == tigrinc::PROMER_STRING)
        this->datatype = tigrinc::PROMER_DATA;
    else
        this->datatype = tigrinc::NULL_DATA;
}

void DeltaFilter::add_record_header(const std::string& idR, const std::string& idQ, long lenR, long lenQ)
{
    if(lenR <= 0 || lenQ <= 0)
    {
        cerr << "lenR <= 0 or lenQ <= 0" << endl;
        exit(1);
    }

    if(!is_first_record)
        insert_dep();
    else
        is_first_record = false;

    record_m.idR = idR;
    record_m.idQ = idQ;
    record_m.lenR = lenR;
    record_m.lenQ = lenQ;

}

void DeltaFilter::add_alignment(long sR, long eR, long sQ, long eQ, long idyc, long simc, long stpc, const vector<long>& all_deltas)
{
    
    if ( sR <= 0  ||  eR <= 0  || sQ <= 0  ||  eQ <= 0  || idyc < 0  ||  simc < 0  ||  stpc < 0 )
    {
        cerr << "sR | eR | sQ | eQ | <=0 or  idyc | simc | stpc < 0" << endl;
        cerr << sR << ' ' << eR << ' ' << sQ << ' ' << eQ << ' ' << idyc << ' ' << simc << ' ' << stpc << endl;
        exit(1);
    }

    record_m.aligns.push_back( tigrinc::DeltaAlignment_t() );
    tigrinc::DeltaAlignment_t& align = record_m.aligns.back();

    align.sR = sR;
    align.eR = eR;
    align.sQ = sQ;
    align.eQ = eQ;
    align.idyc = idyc;
    align.simc = simc;
    align.stpc = stpc;
    align.deltas.reserve( all_deltas.size() + 1 );

    float total = labs(eR - sR) + 1.0;

    if(datatype == tigrinc::PROMER_DATA)
        total /= 3.0;
    
    for(vector<long>::const_iterator d_it = all_deltas.begin(); d_it != all_deltas.end(); ++d_it)
    {
        if(*d_it < 0)
            ++total;
        align.deltas.push_back( *d_it );
    }
    align.deltas.push_back( 0 );

    align.idy = (total - (float)align.idyc) / total * 100.0;
    align.sim = (total - (float)align.simc) / total * 100.0;
    align.stp = (float)align.stpc / (total * 2.0) * 100.0;
}

void DeltaFilter::add_alignment(long sR, long eR, long sQ, long eQ, long idyc, long simc, long stpc)
{
    if ( sR <= 0  ||  eR <= 0  || sQ <= 0  ||  eQ <= 0  || idyc < 0  ||  simc < 0  ||  stpc < 0 )
    {
        cerr << "sR | eR | sQ | eQ | <=0 or  idyc | simc | stpc < 0" << endl;
        cerr << sR << ' ' << eR << ' ' << sQ << ' ' << eQ << ' ' << idyc << ' ' << simc << ' ' << stpc << endl;
        exit(1);
    }

    if(!record_m.aligns.empty())
    {
        tigrinc::DeltaAlignment_t& ali = record_m.aligns.back();
        ali.deltas.push_back( 0 );

        ali.idy = (record_m_total - (float)ali.idyc) / record_m_total * 100.0;
        ali.sim = (record_m_total - (float)ali.simc) / record_m_total * 100.0;
        ali.stp = (float)ali.stpc / (record_m_total * 2.0) * 100.0;
    }

    record_m.aligns.push_back( tigrinc::DeltaAlignment_t() );
    tigrinc::DeltaAlignment_t& align = record_m.aligns.back();

    align.sR = sR;
    align.eR = eR;
    align.sQ = sQ;
    align.eQ = eQ;
    align.idyc = idyc;
    align.simc = simc;
    align.stpc = stpc;

    record_m_total = labs(eR - sR) + 1.0;

    if(datatype == tigrinc::PROMER_DATA)
        record_m_total /= 3.0;
}

void DeltaFilter::add_indel(long d)
{
    tigrinc::DeltaAlignment_t& align = record_m.aligns.back();
    if(d < 0)
        ++record_m_total;
    align.deltas.push_back( d );
}

void DeltaFilter::add_finished()
{
    if(record_m.aligns.empty())
        return;
    tigrinc::DeltaAlignment_t& ali = record_m.aligns.back();
    ali.deltas.push_back( 0 );

    ali.idy = (record_m_total - (float)ali.idyc) / record_m_total * 100.0;
    ali.sim = (record_m_total - (float)ali.simc) / record_m_total * 100.0;
    ali.stp = (float)ali.stpc / (record_m_total * 2.0) * 100.0;
}

void DeltaFilter::run_delta_filter()
{
    if(!is_first_record)
        insert_dep();
    else
        return;

    if ( OPT_MinIdentity > 0  ||  OPT_MinLength > 0 )
        flagScore(OPT_MinLength, OPT_MinIdentity);

    //-- Uniqueness requirements
    if ( OPT_MinUnique > 0 )
        flagUNIQ(OPT_MinUnique);

    //-- Query-based LIS
    if ( OPT_QLIS )
        flagQLIS(OPT_Epsilon, OPT_MaxOverlap);

    //-- Reference-based LIS
    if ( OPT_RLIS )
        flagRLIS(OPT_Epsilon, OPT_MaxOverlap);

    //-- Global LIS
    if ( OPT_GLIS )
        flagGLIS(OPT_Epsilon);

    //-- M-to-M
    if ( OPT_MtoM )
        flagMtoM(OPT_Epsilon, OPT_MaxOverlap);

    //-- 1-to-1
    if ( OPT_1to1 )
        flag1to1(OPT_Epsilon, OPT_MaxOverlap);
}

void DeltaFilter::get_result(delta_filter_ProcessResult& processFilter)
{
    if(is_first_record)
        return;
    bool header;
    long s1, e1, s2, e2;

    map<string, tigrinc::DeltaNode_t>::const_iterator mi;
    vector<tigrinc::DeltaEdge_t *>::const_iterator ei;
    vector<tigrinc::DeltaEdgelet_t *>::const_iterator eli;

    for ( mi = qrynodes.begin(); mi != qrynodes.end(); ++ mi )
    {
        for ( ei  = (mi->second).edges.begin(); ei != (mi->second).edges.end(); ++ ei )
        {
            header = false;

            for ( eli  = (*ei)->edgelets.begin(); eli != (*ei)->edgelets.end(); ++ eli )
            {
                if ( ! (*eli)->isGOOD )
                    continue;

                //-- Print the sequence header
                if ( ! header )
                {
                    processFilter.store_header( processFilter.res,
                            (*ei)->refnode->id->c_str(),
                            (*ei)->qrynode->id->c_str(),
                            (*ei)->refnode->len,
                            (*ei)->qrynode->len);
                    header = true;
                }
                //-- Print the alignment
                s1 = (*eli)->loR;
                e1 = (*eli)->hiR;
                s2 = (*eli)->loQ;
                e2 = (*eli)->hiQ;
                if ( (*eli)->dirR == tigrinc::REVERSE_DIR )
                    swap (s1, e1);
                if ( (*eli)->dirQ == tigrinc::REVERSE_DIR )
                    swap (s2, e2);

                processFilter.new_cluster(processFilter.res, s1, e1, s2, e2,
                        (*eli)->idyc, (*eli)->simc, (*eli)->stpc);

                istringstream iss((*eli)->delta);
                long t_dt;
                while(iss >> t_dt)
                {
                    if(t_dt == 0)
                        break;
                    processFilter.add_delta(processFilter.res, t_dt);
                }
            }
        }
    }
}

void DeltaFilter::clear()
{
    is_first_record = true;
    tigrinc::DeltaGraph_t::clear();
}

