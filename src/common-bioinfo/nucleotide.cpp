#include "nucleotide.h"

namespace loon
{

int Nucleotide::base_to_id(char ch) const
{
    switch(ch)
    {
        case 'A': case 'a':
            return 0;
        case 'T': case 't':
            return 1;
        case 'G': case 'g':
            return 2;
        case 'C': case 'c':
            return 3;
        case '-':
            return 4;
        default:
            ERROR_PRINT("unknown nucleotide(%d)", int(ch));
            exit(1);
    }
}

char Nucleotide::id_to_base(int id) const
{
    static char str[] = "ATGC-";
    if(id < 5 && id >=0)
        return str[id];
    ERROR_PRINT("unknown index %d for base", id);
    exit(1);
}

Nucleotide::Nucleotide(): max_cnt(0), total_cnt(0)
{
    memset(bases_cnt, 0, sizeof(bases_cnt));
}

Nucleotide::Nucleotide(char c, bool count/*=false*/)
{
    memset(bases_cnt, 0, sizeof(bases_cnt));
    init_base(c, count);
}

void Nucleotide::init_base(char c, bool count/*=true*/)
{
    ch = c;
    if(count)
    {
        max_cnt = 1;
        total_cnt = 1;
        bases_cnt[ base_to_id(ch) ] = 1;
    }
    else
    {
        max_cnt = 0;
        total_cnt = 0;
    }
}

void Nucleotide::init_base(char c, int c0, int c1, int c2, int c3, int c4)
{
    ch = c;
    bases_cnt[0] = c0;
    bases_cnt[1] = c1;
    bases_cnt[2] = c2;
    bases_cnt[3] = c3;
    bases_cnt[4] = c4;
    max_cnt = bases_cnt[ base_to_id(c) ];
    total_cnt = c0 + c1 + c2 + c3 + c4;
}

void Nucleotide::init_base(char c, int* cnt)
{
    ch = c;
    memmove(bases_cnt, cnt, sizeof(bases_cnt));
    max_cnt = bases_cnt[ base_to_id(c) ];
    total_cnt = 0;
    for(int i=0; i < 5; ++i)
        total_cnt += cnt[i];
}

char Nucleotide::get_max_base() const
{
    return ch;
}

void Nucleotide::set_base(char c)
{
    ++total_cnt;
    int tmp = (++ bases_cnt[ base_to_id(c) ]);
    if(tmp > max_cnt)
    {
        max_cnt = tmp;
        ch = c;
    }
}

bool Nucleotide::equal(char c) const
{
    return (max_cnt == bases_cnt[ base_to_id(c) ]);
}

bool Nucleotide::simple_equal(char c) const
{
    return (ch == c);
}

bool Nucleotide::equal_with_threshold(char c, double threshold) const
{
    return (max_cnt ? (double(bases_cnt[ base_to_id(c) ])/max_cnt >= threshold) : (c == ch));
}

bool Nucleotide::is_primal() const
{
    return !max_cnt;
}

bool Nucleotide::is_dash() const
{
    return (ch == '-');
}

double Nucleotide::get_ratio(char c) const
{
    return double(bases_cnt[ base_to_id(c) ]) / total_cnt;
}

}
