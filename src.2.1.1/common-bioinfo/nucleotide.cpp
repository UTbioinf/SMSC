#include "nucleotide.h"

namespace loon
{
char Nucleotide::to_upper(char ch) const
{
    if(ch >= 'a' && ch <= 'z')
        return (ch - 'a' + 'A');
    return ch;
}

char* Nucleotide::expand_base(char ch) const
{
    static char str[6] = {0};
    ch = to_upper(ch);
    switch(ch)
    {
        case 'A':
            str[0] = 'A', str[1] = 0;
            break;
        case 'T':
            str[0] = 'T', str[1] = 0;
            break;
        case 'G':
            str[0] = 'G', str[1] = 0;
            break;
        case 'C':
            str[0] = 'C', str[1] = 0;
            break;
        case 'R':
            str[0] = 'A', str[1] = 'G', str[2] = 0;
            break;
        case 'Y':
            str[0] = 'C', str[1] = 'T', str[2] = 0;
            break;
        case 'S':
            str[0] = 'C', str[1] = 'G', str[2] = 0;
            break;
        case 'W':
            str[0] = 'A', str[1] = 'T', str[2] = 0;
            break;
        case 'M':
            str[0] = 'A', str[1] = 'C', str[2] = 0;
            break;
        case 'K':
            str[0] = 'G', str[1] = 'T', str[2] = 0;
            break;
        case 'B':
            str[0] = 'G', str[1] = 'T', str[2] = 'C', str[3] = 0;
            break;
        case 'D':
            str[0] = 'G', str[1] = 'T', str[2] = 'A', str[3] = 0;
            break;
        case 'H':
            str[0] = 'C', str[1] = 'T', str[2] = 'A', str[3] = 0;
            break;
        case 'V':
            str[0] = 'G', str[1] = 'C', str[2] = 'A', str[3] = 0;
            break;
        case 'N':
            str[0] = 'A', str[1] = 'T', str[2] = 'G', str[3] = 'C', str[4] = 0;
            break;
        case '-':
            str[0] = '-', str[1] = 0;
            break;
        default:
            ERROR_PRINT("Unknown base [%d: %c]", int(ch), ch);
            exit(1);
    }
    return str;
}

char Nucleotide::rand_base(char ch) const
{
    ch = to_upper( ch );
    int rand_choice = rand();
    switch(ch)
    {
        case 'A':   return 'A';
        case 'T':   return 'T';
        case 'G':   return 'G';
        case 'C':   return 'C';
        case 'R':   return (( rand_choice & 1 ) ? 'A' : 'G');
        case 'Y':   return (( rand_choice & 1 ) ? 'C' : 'T');
        case 'S':   return (( rand_choice & 1 ) ? 'C' : 'G');
        case 'W':   return (( rand_choice & 1 ) ? 'A' : 'T');
        case 'M':   return (( rand_choice & 1 ) ? 'A' : 'C');
        case 'K':   return (( rand_choice & 1 ) ? 'G' : 'T');
        case 'B':   return id_to_base(rand_choice % 3 + 1);
        case 'D':   return id_to_base(rand_choice % 3);
        case 'H':   return id_to_base((rand_choice % 3) ^ 1);
        case 'V':   return id_to_base((rand_choice % 3) ^ 2);
        case 'N':   return id_to_base(rand_choice & 4);
        case '-':   return '-';
        default:
            ERROR_PRINT("Unknown base [%d: %c]", int(ch), ch);
            exit(1);
    }
}

int Nucleotide::base_to_id(char ch) const
{
    ch = to_upper(ch);
    switch(ch)
    {
        case 'A':
            return 0;
        case 'T':
            return 1;
        case 'G':
            return 2;
        case 'C':
            return 3;
        case '-':
            return 4;
        default:
            ERROR_PRINT("unknown nucleotide[%d, %c]", int(ch), ch);
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
    srand(time(NULL));
    memset(bases_cnt, 0, sizeof(bases_cnt));
}

Nucleotide::Nucleotide(char c, bool count/*=false*/)
{
    srand(time(NULL));
    memset(bases_cnt, 0, sizeof(bases_cnt));
    init_base(c, count);
}

void Nucleotide::init_base(char c, bool count/*=true*/)
{
    ch = rand_base(c);
    if(count)
    {
        max_cnt = 1;
        total_cnt = 0;
        const char* str = expand_base( ch );
        for(; str[ total_cnt ]; ++total_cnt)
        {
            bases_cnt[ base_to_id( str[total_cnt] ) ] = 1;
        }
    }
    else
    {
        max_cnt = 0;
        total_cnt = 0;
    }
}

void Nucleotide::init_base(char c, int c0, int c1, int c2, int c3, int c4)
{
    ch = to_upper(c);
    bases_cnt[0] = c0;
    bases_cnt[1] = c1;
    bases_cnt[2] = c2;
    bases_cnt[3] = c3;
    bases_cnt[4] = c4;
    max_cnt = bases_cnt[ base_to_id(ch) ];
    total_cnt = c0 + c1 + c2 + c3 + c4;
}

void Nucleotide::init_base(char c, int* cnt)
{
    ch = to_upper(c);
    memmove(bases_cnt, cnt, sizeof(bases_cnt));
    max_cnt = bases_cnt[ base_to_id(ch) ];
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
    const char* str = expand_base( ch );
    int old_ch = ch;
    for(int i = 0; str[i]; ++i)
    {
        ++total_cnt;
        int tmp = (++ bases_cnt[ base_to_id( str[i] ) ]);
        if(tmp > max_cnt || (tmp == max_cnt && str[i] == old_ch))
        {
            max_cnt = tmp;
            ch = str[i];
        }
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
    if(total_cnt == 0)
        return 0;
    return double(bases_cnt[ base_to_id(c) ]) / total_cnt;
}

}
