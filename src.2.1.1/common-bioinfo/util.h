#ifndef LOON_UTIL_H
#define LOON_UTIL_H

#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>

#include "baseheader.h"

namespace loon
{

void open_file(std::ifstream& fin, const std::string& fname);
void open_file(std::ofstream& fout, const std::string& fname);
void check_file_open(std::ifstream& fin, const std::string& fname);
void check_file_open(std::ofstream& fout, const std::string& fname);
void compute_reverse_complement(const std::string& str, std::string& dest_str);
void compute_reverse_complement(std::string& str);
void save_fasta(const std::string& tag, const std::string& seq, const std::string& fname);

}
#endif
