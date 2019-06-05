#ifndef __MSA_CONSENSUS_H
#define __MSA_CONSENSUS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>

namespace loon
{

// all locations are left included, right excluded
class PairwiseAln
{
public:
    size_t ref_start, ref_end;
    size_t qry_start, qry_end;
    std::string* qry_str;
    std::vector<long long> cigar_cnt;
    std::string cigar_mttype;
    size_t cigar_loc, ref_loc, qry_loc;
public:
    void clear();
    int next_cigar();
    char get_cigar_type() const;
    long long get_cigar_cnt() const;
    std::string get_deletion_str() const;
    char get_qry_char() const;
    void add_alignment(char& cur_state, char match_type, long long cnt = 1);

    bool operator<(const PairwiseAln& rhs) const;
};

class MSA_Consensus
{
public:
    const static char* nucs;
    std::string* ref_str;
    std::string* ref_quality;
    std::vector<PairwiseAln> alns;
    std::vector<std::string> cons;
    std::vector<std::string> cons_score;
public:
    void add_ref(std::string* ref_str, std::string* ref_quality);
    void clear();
    void reserve(size_t n);
    void add_alignment(size_t ref_start, size_t ref_end, size_t qry_start, size_t qry_end,
            std::string& qry_str,
            std::vector<long long>& cigar_cnt, std::string& cigar_mttype); // cigar_cnt and cigar_mttype will be swapped
    void add_alignment(size_t ref_start, size_t ref_end, size_t qry_start, size_t qry_end,
            std::string& qry_str,
            std::vector<long>& indels);
    void sort();
    void do_consensus(bool include_ref = false);
    std::vector<std::string>& get_consensus();
    std::vector<std::string>& get_consensus_score();
    char accuracyToAscii(double score);

    void remove_D_ids(std::vector<size_t>& MI_ids, std::vector<size_t>& D_ids);
    void consensus_deletion(const std::vector<std::string>& msa_aln, size_t total);
    void consensus_deletion(std::vector< std::vector<size_t> >& nuc_cnt, size_t total, size_t extra_gap);
    size_t max_id(long long a[], size_t n);
    size_t max_id(std::vector<size_t>& a);
    void inc_nuc_cnt(long long nuc_cnt[], char nuc_type);

    void write_for_debug();
    void read_for_debug();
};

} // end of namespace loon
#endif
