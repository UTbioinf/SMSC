#ifndef _MUMMER_API_H
#define _MUMMER_API_H

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <errno.h>
#include <unistd.h>
#include <cstring>
#include <limits.h>


#include "nucmer.h"
#include "common-bio/Consensus.h"
#include "common-bio/Multiseq.h"
#include "common-bio/util.h"
#include "common-bio/baseheader.h"
#include "common-bio/Progress.h"

#ifdef __APPLE__
    #include <libproc.h>
#endif

//#define DONTCLEAN

namespace loon
{
#define MUMMER_ENGINE 0
#define BLASR_ENGINE 1

class MUMTaskParams
{
public:
    const std::string& ref_file;
    const std::string& qry_file;
    const std::string& rev_qry_file;
    const std::string prefix;
    const std::vector<size_t>& ref_starts;
    const Multiseq* refp;
    const Multiseq* qryp;
    double diagfactor;
    int min_match_length, min_cluster, max_gap, diagdiff, breaklen;
    bool remove_middle;
public:
    MUMTaskParams(const std::string& rf, const std::string& qf,
            const std::string& rev_qf, const std::string& pfx,
            const std::vector<size_t>& ref_s = std::vector<size_t>(),
            bool rmd = false,
            const Multiseq* rp = NULL, const Multiseq* qp = NULL,
            int min_ml = 20, int min_c = 200,
            int max_g = 20, int ddiff = 5,
            double diagf = .12, int brklen = 150):
                ref_file(rf), qry_file(qf),
                rev_qry_file(rev_qf), prefix(pfx),
                ref_starts( ref_s ),
                refp(rp), qryp(qp),
                diagfactor(diagf),
                min_match_length( min_ml), min_cluster(min_c),
                max_gap(max_g), diagdiff(ddiff),
                breaklen(brklen),
                remove_middle(rmd)
    {}
};

void mumapi_set_script_path();
void mumapi_unsetdir();
void mumapi_setdir(const std::string& path);
void mumapi_mkdir();
void mumapi_rmdir();
void mumapi_clean_exit(int status);
void mumapi_setpoolsize(int ps);

void mumapi_set_engine(int engine);

template<class T>
char* const to_charp(T num, std::string& str)
{
    std::ostringstream oss;
    oss << num;
    str = oss.str();
    return (char *const)str.c_str();
}

void make_ref_intervals(const Multiseq& ref, std::vector<size_t>& ref_starts);
size_t get_interval(const std::vector<size_t>& ref_starts, size_t pos);
bool keepthis(const std::vector<size_t>& ref_starts, const std::string& line, 
        size_t qlen);
void remove_middle_matching(const std::vector<size_t>& ref_starts,
        const Multiseq& qry,
        const std::string& infile, const std::string& outfile);

void exec_external_command_with_pipe(const char* path1, char *const args1[],
        const char* path2, char *const args2[]);
void exec_external_command(const char* path, char *const args[]);
void run_nucmer_task(void* arg);
        
void read_deltas(nucmer_Deltas& deltas, const std::string& fname);
void all_v_all_matching(const Multiseq& ref, const Multiseq& qry,
        nucmer_Deltas& deltas_forward, nucmer_Deltas& deltas_backward,
        bool ends_matching= false);

void all_v_all_matching_mummer(const Multiseq& ref, const Multiseq& qry, 
        nucmer_Deltas& deltas_forward, nucmer_Deltas& deltas_backward,
        bool ends_matching = false);

// api for blasr
void all_v_all_matching_blasr(const Multiseq& ref, const Multiseq& qry, 
        nucmer_Deltas& deltas_forward, nucmer_Deltas& deltas_backward,
        bool ends_matching = false);

void blasr_align_init(const std::string& X, const std::string& Y);
void blasr_run();
void blasr_output_transform_forward();

/*
choice:
    * 3: align_full
    * 1: align_right
    * 2: align_left
*/
void blasr_align_choice(const std::string& X, const std::string& Y,
        nucmer_Deltas& deltas, size_t boundary, int choice);
}

#endif
