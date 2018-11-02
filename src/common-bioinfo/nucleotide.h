#ifndef __NUCLEOTIDE_H
#define __NUCLEOTIDE_H

#include <iostream>
#include <vector>

#include <cstring>
#include <cstdlib>

#include "baseheader.h"

namespace loon
{

class Nucleotide
{
public:
    char ch;
    int max_cnt, total_cnt;
    int bases_cnt[5];

    int base_to_id(char ch) const;
    char id_to_base(int id) const;
public:
    Nucleotide();
    Nucleotide(char c, bool count=false);
    void init_base(char c, bool count = false); // if count == true: bases_cnt[c] = 1
    void init_base(char c, int c0, int c1, int c2, int c3, int c4);
    void init_base(char c, int* cnt);
    char get_max_base() const;
    void set_base(char c); // W/ ++bases_cnt[c]
    bool equal(char c) const;
    bool simple_equal(char c) const;
    bool equal_with_threshold(char c, double threshold) const;
    bool is_primal() const;
    bool is_dash() const;
    double get_ratio(char c) const;
};

}

#endif
