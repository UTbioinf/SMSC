#ifndef LOON_MULTISEQ_H
#define LOON_MULTISEQ_H

#include<vector>
#include<string>

namespace loon
{

class Fasta
{
public:
    std::string tag;
    std::string tag_remaining;
    std::string seq;
    std::string rev;
public:
    char base_rc(char c);
    void save_fasta(std::ostream& out);
    void compute_reverse_complement();
};

class Multiseq: public std::vector<Fasta>
{
public:
    void read_fasta(const char* fname);
    void save_fasta(const char* fname);
    void compute_reverse_complement();
};

}

#endif
