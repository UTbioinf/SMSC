#ifndef __CONSENSUS_H
#define __CONSENSUS_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <cstdlib>

#include "baseheader.h"
#include "nucleotide.h"
#include "nucmer.h"
#include "util.h"
#include "Aligner.h"

namespace loon
{// namespace loon starts

class Consensus: public std::vector<Nucleotide>
{
public:
    std::string str;
    size_t left_delimiter;
    size_t right_delimiter;
    static Aligner aligner;

public:
    Consensus();
    void swap(Consensus& t);
    void swap(std::vector<Nucleotide>& t);
    void set_str(const std::string& str);
    std::string& toString();
    void toString(std::string& str);
    void consensus_init(const std::string& qry);
    void consensus_error(const std::string& qry);
    void consensus_result(const std::string& qry);

    void consensus_full(const std::string& str, bool unaligned_stop = false);
    void consensus_full(void (*compute_deltas_func)(const std::string&,
                    const std::string&, nucmer_Deltas&, size_t, int), 
            const std::string& str, bool unaligned_stop = false);

    void consensus_left(const std::string& str, size_t right_most = INF, 
            bool unaligned_stop = false);
    void consensus_left(void (*compute_deltas_func)(const std::string&,
                    const std::string&, nucmer_Deltas&, size_t, int),
            const std::string& str, 
            size_t right_most = INF, bool unaligned_stop = false);

    void consensus_right(const std::string& str, size_t left_most = INF,
            bool unaligned_stop = false);
    void consensus_right(void (*compute_deltas_func)(const std::string&,
                    const std::string&, nucmer_Deltas&, size_t, int),
            const std::string& str,
            size_t left_most = INF, bool unaligned_stop = false);

    void consensus_left_auto_delimiter(const std::string& str, 
            bool unaligned_stop = false);
    void consensus_left_auto_delimiter(void (*compute_deltas_func)(const std::string&,
                    const std::string&, nucmer_Deltas&, size_t, int),
            const std::string& str, bool unaligned_stop = false);

    void consensus_right_auto_delimiter(const std::string& str, 
            bool unaligned_stop = false);
    void consensus_right_auto_delimiter(void (*compute_deltas_func)(const std::string&,
                    const std::string&, nucmer_Deltas&, size_t, int),
            const std::string& str, bool unaligned_stop = false);

    void save_binary(std::ostream& out);
    void save_binary(const std::string& fname);
    void read_binary(std::istream& in);
    void read_binary(const std::string& fname);
};

}// namespace loon ends
#endif
