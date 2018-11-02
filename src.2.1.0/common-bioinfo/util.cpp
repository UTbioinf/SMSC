#include <cstdlib>
#include "util.h"

namespace loon
{

void open_file(std::ifstream& fin, const std::string& fname)
{
    fin.open(fname.c_str());
    check_file_open(fin, fname);
}

void open_file(std::ofstream& fout, const std::string& fname)
{
    fout.open(fname.c_str());
    check_file_open(fout, fname);
}

void check_file_open(std::ifstream& fin, const std::string& fname)
{
    if(!fin.is_open())
    {
        ERROR_PRINT("Cannot open file \"%s\"", fname.c_str());
        exit(1);
    }
}

void check_file_open(std::ofstream& fout, const std::string& fname)
{
    if(!fout.is_open())
    {
        ERROR_PRINT("Cannot open file \"%s\"", fname.c_str());
        exit(1);
    }
}


static char base_rc(char c)
{
    if(c >= 'a' && c <= 'z')
        c = c - 'a' + 'A';
    switch(c)
    {
        case 'A':   return 'T';
        case 'T':   return 'A';
        case 'G':   return 'C';
        case 'C':   return 'G';
        case 'R':   return 'Y';
        case 'Y':   return 'R';
        case 'S':   return 'S';
        case 'W':   return 'W';
        case 'M':   return 'K';
        case 'K':   return 'M';
        case 'B':   return 'V';
        case 'D':   return 'H';
        case 'H':   return 'D';
        case 'V':   return 'B';
        case 'N':   return 'N';
        case '-':   return '-';
        default:
            ERROR_PRINT("Unknown base [%d: %c]", int(c), c);
            exit(1);
    };
}

void compute_reverse_complement(const std::string& str, std::string& dest_str)
{
    dest_str = str;

    char c;
    size_t i, j;
    for(i = 0, j = dest_str.length() - 1; i < j; ++i, --j)
    {
        c = base_rc( dest_str[i] );
        dest_str[i] = base_rc( dest_str[j] );
        dest_str[j] = c;
    }

    if(i == j)
        dest_str[i] = base_rc( dest_str[i] );
}

void compute_reverse_complement(std::string& str)
{
    char c;
    size_t i, j;
    for(i = 0, j = str.length() - 1; i < j; ++i, --j)
    {
        c = base_rc( str[i] );
        str[i] = base_rc( str[j] );
        str[j] = c;
    }

    if(i == j)
        str[i] = base_rc( str[i] );
}

void save_fasta(const std::string& tag, const std::string& seq, const std::string& fname)
{
    std::ofstream fout(fname.c_str());
    check_file_open(fout, fname);

    fout << ">" << tag << std::endl;
    for(size_t i = 0; i < seq.length(); i += 80)
        fout << seq.substr(i, 80) << std::endl;
    fout.close();
}
}
