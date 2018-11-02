#include <fstream>
#include <cctype>
#include "Multiseq.h"
#include "util.h"

namespace loon
{

void Fasta::save_fasta(std::ostream& out)
{
    out << ">" << tag << tag_remaining << std::endl;
    for(size_t i = 0; i < seq.length(); i += 80)
        out << seq.substr(i, 80) << std::endl;
}

char Fasta::base_rc(char c)
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
    }
}

void Fasta::compute_reverse_complement()
{
    rev.clear();
    rev.reserve( seq.length() );
    for(std::string::reverse_iterator rit = seq.rbegin(); rit != seq.rend(); ++rit)
        rev.push_back( base_rc( *rit ) );
}

void Multiseq::read_fasta(const char* fname)
{
    std::ifstream fin(fname);
    check_file_open(fin, fname);

    std::string line;
    while(std::getline(fin, line))
    {
        if(line[0] == '>')
        {
            this->push_back( Fasta() );
            Fasta& fasta = this->back();

            std::size_t space_pos = line.find(' ');
            if(space_pos != std::string::npos)
            {
                fasta.tag = line.substr(1, space_pos-1);
                fasta.tag_remaining = line.substr(space_pos);
            }
            else
                fasta.tag = line.substr(1);
        }
        else if(line[0] != ';')
        {
            for(std::size_t i = 0; i<line.length(); ++i)
                line[i] = toupper(line[i]);
            this->back().seq += line;
        }
    }
    fin.close();
}

void Multiseq::save_fasta(const char* fname)
{
    std::ofstream fout(fname);
    check_file_open(fout, fname);

    for(std::vector<Fasta>::iterator it = this->begin(); it != this->end(); ++it)
    {
        fout << ">" << (it->tag) << (it->tag_remaining) << std::endl;
        for(size_t i = 0; i < it->seq.length(); i += 80)
            fout << (it->seq.substr(i, 80)) << std::endl;
    }

    fout.close();
}

void Multiseq::compute_reverse_complement()
{
    for(std::vector<Fasta>::iterator it = this->begin(); it != this->end(); ++it)
        it->compute_reverse_complement();
}

}
