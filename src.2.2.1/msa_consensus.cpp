#include <cmath>
#include "msa_consensus.h"
#include "mummer_api.h"

namespace loon
{
void PairwiseAln::clear()
{
    cigar_cnt.clear();
    cigar_mttype.clear();
}

// return value:
// * -1: finished processing this alignment
// *  0: cigar_mttype not change
// *  1: cigar_mttype changed
int PairwiseAln::next_cigar()
{
    if(cigar_mttype[ cigar_loc ] == 'D')
    {
        qry_loc += cigar_cnt[ cigar_loc ];
        ++cigar_loc;
        if(cigar_loc == cigar_mttype.length())  return -1;
        else    return 1;
    }
    else
    {
        if(cigar_mttype[ cigar_loc ] == 'I')
        {
            ++ref_loc;
            if( --cigar_cnt[ cigar_loc ] == 0)
            {
                ++cigar_loc;
                if(cigar_loc == cigar_cnt.size())   return -1;
                return 1;
            }
            return 0;
        }
        else // M, X or =
        {
            ++ref_loc;  ++qry_loc;
            if( --cigar_cnt[ cigar_loc ] == 0)
            {
                ++cigar_loc;
                if(cigar_loc == cigar_cnt.size())   return -1;
                return 1;
            }
            return 0;
        }
    }
}

char PairwiseAln::get_cigar_type() const
{
    return cigar_mttype[ cigar_loc ];
}

long long PairwiseAln::get_cigar_cnt() const
{
    return cigar_cnt[ cigar_loc ];
}

std::string PairwiseAln::get_deletion_str() const
{
    return qry_str->substr( qry_loc, cigar_cnt[ cigar_loc ] );

}

char PairwiseAln::get_qry_char() const
{
    if(cigar_mttype[ cigar_loc ] == 'I')
        return '-';
    return qry_str->at( qry_loc );
}

void PairwiseAln::add_alignment(char& cur_state, char match_type, long long cnt /* = 1*/)
{
    if(cur_state == match_type)
        cigar_cnt.back() += cnt;
    else
    {
        cigar_mttype.push_back( match_type );
        cigar_cnt.push_back( cnt );
        cur_state = match_type;
    }
}

bool PairwiseAln::operator<(const PairwiseAln& rhs) const
{
    return ref_start < rhs.ref_start;
}




const char* MSA_Consensus::nucs = "ATGC-";

void MSA_Consensus::add_ref(std::string* ref_str, std::string* ref_quality)
{
    this->ref_str = ref_str;
    this->ref_quality = ref_quality;
}

void MSA_Consensus::clear()
{
    alns.clear();
    cons.clear();
    cons_score.clear();
}

void MSA_Consensus::reserve(size_t n)
{
    alns.reserve( n );
}

void MSA_Consensus::add_alignment(size_t ref_start, size_t ref_end, size_t qry_start, size_t qry_end,
        std::string& qry_str,
        std::vector<long long>& cigar_cnt, std::string& cigar_mttype)
{
    alns.push_back( PairwiseAln() );
    PairwiseAln& tmp_aln = alns.back();
    tmp_aln.ref_start = ref_start;
    tmp_aln.ref_end = ref_end;
    tmp_aln.qry_start = qry_start;
    tmp_aln.qry_end = qry_end;
    tmp_aln.qry_str = &qry_str;
    tmp_aln.cigar_cnt.swap( cigar_cnt );
    tmp_aln.cigar_mttype.swap( cigar_mttype );
    tmp_aln.cigar_loc = 0;
    tmp_aln.ref_loc = ref_start;
    tmp_aln.qry_loc = qry_start;
}

void MSA_Consensus::add_alignment(size_t ref_start, size_t ref_end, size_t qry_start, size_t qry_end,
        std::string& qry_str, std::vector<long>& indels)
{
    alns.push_back( PairwiseAln() );
    PairwiseAln& tmp_aln = alns.back();
    tmp_aln.ref_start = ref_start;
    tmp_aln.ref_end = ref_end;
    tmp_aln.qry_start = qry_start;
    tmp_aln.qry_end = qry_end;
    tmp_aln.qry_str = &qry_str;
    tmp_aln.cigar_loc = 0;
    tmp_aln.ref_loc = ref_start;
    tmp_aln.qry_loc = qry_start;

    tmp_aln.cigar_cnt.clear();
    tmp_aln.cigar_mttype.clear();
    char cur_state = 0;
    size_t pos_x = ref_start, pos_y = qry_start;
    for(std::vector<long>::const_iterator indel_it = indels.begin(); indel_it != indels.end(); ++indel_it)
    {
        long steps = std::abs(*indel_it);
        if(steps > 1)
        {
            tmp_aln.add_alignment(cur_state, 'M', steps - 1);
            pos_x += steps - 1;
            pos_y += steps - 1;
        }
        if(*indel_it > 0)
        {
            tmp_aln.add_alignment(cur_state, 'I');
            ++pos_x;
        }
        else
        {
            tmp_aln.add_alignment(cur_state, 'D');
            ++pos_y;
        }
    }
    if(pos_x < ref_end)
    {
        tmp_aln.add_alignment(cur_state, 'M', ref_end - pos_x);
    }
}

void MSA_Consensus::sort()
{
    std::sort(alns.begin(), alns.end());
}

void MSA_Consensus::do_consensus(bool include_ref/* = false */)
{
    if(alns.empty())    return;
    std::vector<size_t> MI_ids, D_ids;
    size_t cur_aln = 0, cur_ref_loc = alns[ 0 ].ref_start;
    size_t n_alns = alns.size();

    cons.push_back( std::string() );
    cons_score.push_back( std::string() );

    if(include_ref)
    {// deal with unaligned left padding of ref
        if(cur_ref_loc > 0)
        {
            cons.back().assign( ref_str->begin(), ref_str->begin() + cur_ref_loc );
            if(ref_quality) cons_score.back().assign( ref_quality->begin(), ref_quality->begin() + cur_ref_loc );
            else    cons_score.back().assign( cur_ref_loc, '!' );
        }
    }

    std::string muscle_infile = muscle_get_infile_name();
    // deal with the middle aligned part
    while( ! (cur_aln >= n_alns && MI_ids.empty() && D_ids.empty() ))
    {
        // add all the appropriate alignments
        for(; cur_aln < n_alns && alns[ cur_aln ].ref_start <= cur_ref_loc; ++cur_aln)
        {
            char mttype = alns[ cur_aln ].cigar_mttype[ 0 ];
            if( mttype == 'M' || mttype == 'X' || mttype == '=' || mttype == 'I')
                MI_ids.push_back( cur_aln );
            else if(mttype == 'D')
                D_ids.push_back( cur_aln );
            else
            {
                std::cerr << "[" << __FILE__ << ' ' << __LINE__ << "] [ERROR]: "
                        <<  "This implementation only supports M, X, =, I, D in cigar, encountered [" 
                        << mttype << ", " << int(mttype) << "]" << std::endl;
                exit(1);
            }
        }

        // Deal with deletion
        if(!D_ids.empty())
        {
            size_t D_size = D_ids.size();
            size_t total = MI_ids.size() + D_size + (include_ref ? 1 : 0);
            
            if(D_size >= 2 && D_size > (total >> 1))
            {
                std::ofstream fout( muscle_infile.c_str() );
                if(!fout.is_open())
                {
                    std::cerr << "[" << __FILE__ << " " << __LINE__ << "] [ERROR]: Cannot open muscle infile" << std::endl;
                    exit(1);
                }
                for(size_t i = 0; i < D_ids.size(); ++i)
                {
                    fout << ">id_" << i << std::endl;
                    fout << alns[ D_ids[i] ].get_deletion_str() << std::endl;
                }
                fout.close();

                std::vector<std::vector<size_t> > nuc_cnt;
                muscle_compute_msa( nuc_cnt );

                consensus_deletion(nuc_cnt, total, total - D_size);
            }
            remove_D_ids(MI_ids, D_ids);
        }

        if(MI_ids.empty())
        {
            if(cur_aln < n_alns)
            {
                if(include_ref)
                {
                    cons.back().append( ref_str->begin() + cur_ref_loc, ref_str->begin() + alns[ cur_aln ].ref_start );
                    if(ref_quality) cons_score.back().append( ref_quality->begin() + cur_ref_loc, ref_quality->begin() + alns[ cur_aln ].ref_start );
                    else    cons_score.back().append( alns[ cur_aln ].ref_start - cur_ref_loc, '!' );
                }
                else
                {
                    cons.push_back( std::string() );
                    cons_score.push_back( std::string() );
                }
                cur_ref_loc = alns[ cur_aln ].ref_start;
            }
            else    return;
        }
        else
        {
            // Deal with insertion, match and mismatch
            long long nuc_cnt[5] = {0};
            if(include_ref)     inc_nuc_cnt( nuc_cnt, ref_str->at( cur_ref_loc ) );
            for(std::vector<size_t>::iterator m_it = MI_ids.begin(); m_it != MI_ids.end(); ++m_it)
                inc_nuc_cnt( nuc_cnt, alns[ *m_it ].get_qry_char() );
            size_t m_id = max_id(nuc_cnt, 5);
            if(m_id < 4)
            {
                cons.back().push_back( nucs[ m_id ] );
                cons_score.back().push_back( accuracyToAscii( nuc_cnt[m_id] / double( MI_ids.size() + (include_ref ? 1 : 0)) ) );
            }
            size_t n_mi = MI_ids.size();
            for(size_t i = 0; i < n_mi; ++i)
            {
                int ret = alns[ MI_ids[ i ] ].next_cigar();
                if(ret == 0)    continue;
                else if(ret == -1)
                {
                    MI_ids[i] = MI_ids.back();
                    MI_ids.pop_back();
                    --n_mi; --i;
                }
                else
                {
                    if( alns[ MI_ids[i] ].get_cigar_type() == 'D')
                    {
                        D_ids.push_back( MI_ids[i] );
                        MI_ids[ i ] = MI_ids.back();
                        MI_ids.pop_back();
                        --n_mi; --i;
                    }
                }
            }
            ++cur_ref_loc;
        }
    }

    if(include_ref)
    {//deal with unaligned right padding of ref
        if(cur_ref_loc < ref_str->length())
        {
            cons.back().append( ref_str->begin() + cur_ref_loc, ref_str->end());
            if(ref_quality) cons_score.back().append( ref_quality->begin() + cur_ref_loc, ref_quality->end() );
            else    cons_score.back().append(ref_str->length() - cur_ref_loc, '!');
        }
    }
}

std::vector<std::string>& MSA_Consensus::get_consensus()
{
    return cons;
}

std::vector<std::string>& MSA_Consensus::get_consensus_score()
{
    return cons_score;
}

char MSA_Consensus::accuracyToAscii(double score)
{
    int Q = round(-10 * log(1 - score)) + 33;
    if(Q < 33)  return '!';
    if(Q > 126) return '~';
    return char(Q);
}

void MSA_Consensus::remove_D_ids(std::vector<size_t>& MI_ids, std::vector<size_t>& D_ids)
{
    for(std::vector<size_t>::iterator d_it = D_ids.begin(); d_it != D_ids.end(); ++d_it)
        if(alns[ *d_it ].next_cigar() == 1)
        {
            char cigar_type = alns[ *d_it ].get_cigar_type();
            MI_ids.push_back( *d_it );
        }
    D_ids.clear();
}

void MSA_Consensus::consensus_deletion(const std::vector<std::string>& msa_aln, size_t total)
{
    size_t n_aln = msa_aln.size();
    size_t aln_len = msa_aln.front().length();
    for(size_t l = 0; l < aln_len; ++l)
    {
        long long nuc_cnt[5] = {0};
        for(size_t i = 0; i < n_aln; ++i)
            inc_nuc_cnt(nuc_cnt, msa_aln[i][l]);
        nuc_cnt[4] += total - n_aln;
        size_t m_id = max_id(nuc_cnt, 5);
        if(m_id < 4)
        {
            cons.back().push_back( nucs[ m_id ] );
            cons_score.back().push_back( accuracyToAscii( nuc_cnt[ m_id ] / double(total) ) );
        }
    }
}

void MSA_Consensus::consensus_deletion(std::vector<std::vector<size_t> >& nuc_cnt, size_t total, size_t extra_gap)
{
    for(std::vector<std::vector<size_t> >::iterator it = nuc_cnt.begin();
            it != nuc_cnt.end(); ++it)
    {
        it->at(4) += extra_gap;
        size_t m_id = max_id( *it );
        if(m_id < 4)
        {
            cons.back().push_back( nucs[m_id] );
            cons_score.back().push_back( accuracyToAscii( it->at( m_id ) / double(total) ) );
        }
    }
}

size_t MSA_Consensus::max_id(long long a[], size_t n)
{
    size_t max_id = 0;
    for(size_t id = 1; id < n; ++id)
        if(a[ max_id ] < a[ id ])
            max_id = id;
    return max_id;
}

size_t MSA_Consensus::max_id(std::vector<size_t>& a)
{
    size_t max_id = 0;
    for(size_t id = 1; id < a.size(); ++id)
        if(a[ max_id] < a[id] )
            max_id = id;
    return max_id;
}

void MSA_Consensus::inc_nuc_cnt(long long nuc_cnt[], char nuc_type)
{
    switch(nuc_type)
    {
        case 'A':
            ++nuc_cnt[0];
            return;
        case 'T':
            ++nuc_cnt[1];
            return;
        case 'G':
            ++nuc_cnt[2];
            return;
        case 'C':
            ++nuc_cnt[3];
            return;
        case '-':
            ++nuc_cnt[4];
            return;
    }
}

void MSA_Consensus::write_for_debug()
{
    std::ofstream fout("msa_consensus.debug");
    if(ref_quality == NULL)
        fout << 1 << ' ' << (*ref_str) << std::endl;
    else
        fout << 2 << ' ' << (*ref_str) << ' ' << (*ref_quality) << std::endl;
    fout << alns.size() << std::endl;
    for(std::vector<PairwiseAln>::const_iterator alns_it = alns.begin();
            alns_it != alns.end(); ++alns_it)
    {
        size_t nn = alns_it->cigar_cnt.size();
        fout << (*(alns_it->qry_str)) << std::endl;
        fout << (alns_it->ref_start) << ' ' << (alns_it->ref_end) << ' '
                << (alns_it->qry_start) << ' ' << (alns_it->qry_end) << ' '
                << nn << std::endl;
        for(size_t ii = 0; ii < nn; ++ii)
            fout << (alns_it->cigar_cnt[ ii ]) << ' ' << (alns_it->cigar_mttype[ ii ]) << ' ';
        fout << std::endl;
    }
    fout.close();
}

void MSA_Consensus::read_for_debug()
{
    std::ifstream fin("msa_consensus.debug");
    short has_quality;
    fin >> has_quality;
    if(has_quality == 1)
    {
        ref_str = new std::string();
        fin >> (*ref_str);
        ref_quality = NULL;
    }
    else
    {
        ref_str = new std::string(); ref_quality = new std::string();
        fin >> (*ref_str) >> (*ref_quality);
    }
    size_t alns_size;
    fin >> alns_size;
    alns.assign(alns_size, PairwiseAln());

    for(size_t alns_i = 0; alns_i < alns_size; ++alns_i)
    {
        alns[ alns_i ].qry_str = new std::string();
        fin >> *(alns[ alns_i ].qry_str);
        size_t nn;
        fin >> alns[ alns_i ].ref_start >> alns[ alns_i ].ref_end
                >> alns[ alns_i ].qry_start >> alns[ alns_i ].qry_end
                >> nn;
        alns[ alns_i ].cigar_cnt.assign( nn, 0 );
        alns[ alns_i ].cigar_mttype.assign( nn, 0 );
        for(size_t ii = 0; ii < nn; ++ii)
            fin >> alns[ alns_i ].cigar_cnt[ ii ] >> alns[ alns_i ].cigar_mttype[ ii ]; 
                
    }
    fin.close();
}
}// end of namespace loon
