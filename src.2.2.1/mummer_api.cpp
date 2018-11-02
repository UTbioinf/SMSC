#include <sstream>
#include <cmath>
#include <queue>
#include <algorithm>

#include "mummer_api.h"
#include "threadpool.h"

namespace loon
{

template<class T>
char* const to_charp(T num, std::string& str)
{
    std::ostringstream oss;
    oss << num;
    str = oss.str();
    return (char *const)str.c_str();
}

class FileSizeElement
{
public:
    unsigned long long file_size;
    size_t file_idx;
public:
    FileSizeElement(unsigned long long fs = 0, size_t fi = INF):
        file_size(fs), file_idx(INF)
    {}
    void add_file_size(unsigned long long fs)
    {
        file_size += fs;
    }
    void set_file_idx(size_t fi)
    {
        file_idx = fi;
    }
    bool operator<(const FileSizeElement& f) const
    {
        return (file_size > f.file_size || (file_size == f.file_size && file_idx < f.file_idx));
    }
};

static std::string build_default_path()
{
    std::ostringstream oss;
    oss << "/tmp/mummer_api_" << getpid();
    return oss.str();
}

static string get_bin_path()
{
#ifdef __APPLE__
    pid_t pid = getpid();
    char exec_path[ PROC_PIDPATHINFO_MAXSIZE ];
    int ret = proc_pidpath(pid, exec_path, PROC_PIDPATHINFO_MAXSIZE);
    if(ret <= 0)
    {
        fprintf(stderr, "PID %d: proc_pidpath();\n", pid);
        fprintf(stderr, "    %s\n", strerror(errno));
        exit(1);
    }
#elif __linux__
    char exec_path[ PATH_MAX + 1];
    ssize_t ret = readlink("/proc/self/exe", exec_path, PATH_MAX);
    if(ret == -1)
    {
        fprintf(stderr, "readlink error: %s\n", strerror(errno));
        exit(1);
    }
    exec_path[ret] = 0;
#else
    cout << "Unknown operating system, please add the way to get executable path by yourself" << endl;
#endif

    char pathbuf[ PATH_MAX + 1];
    realpath(exec_path, pathbuf);
    char* pch = strrchr(pathbuf, '/');
    pathbuf[ pch - pathbuf ] = 0;
    return string(pathbuf);
}

static std::string tmp_path = build_default_path();
static std::string script_path;
static int poolsize = 1;
static int choose_engine = MUMMER_ENGINE;


static void show_command(const std::string& command)
{
#ifdef DEBUG
    std::cerr << command << std::endl;
#endif
}

static size_t str2lu(const std::string& str)
{
    size_t ret;
    istringstream iss(str);
    iss >> ret;
    return ret;
}

static std::string get_fname(const std::string& prefix, size_t id, const std::string& suffix)
{
    std::string ret;
    ostringstream oss;
    oss << prefix << id << suffix;
    ret = oss.str();
    return ret;
}

static void save_fasta(std::ostream& out, int id, const std::string& seq)
{
    if(seq.empty()) return;
    out << ">" << id << std::endl;
    for(size_t i = 0; i < seq.length(); i += 80)
        out << seq.substr(i, 80) << std::endl;
}

static unsigned long long get_total_length(const Multiseq& multiseq)
{
    unsigned long long ret = 0;
    size_t n = multiseq.size();
    for(size_t i = 0; i < n; ++i)
        ret += multiseq[ i ].seq.length();
    return ret;
}

static void save_ref_files(const Multiseq& ref, std::vector<std::string>& ref_files,
        std::vector<std::vector<size_t> >& ref_start_vector, size_t n, bool ends_matching)
{
    std::ofstream* out_files;
    out_files = new std::ofstream[n];
    std::vector<FileSizeElement> init_queue(n);
    for(size_t i = 0; i < n; ++i)
    {
        // generate file name
        ostringstream oss;
        oss << tmp_path << "/ref." << i << ".fa";
        ref_files[ i ] = oss.str();

        // open file
        out_files[i].open( ref_files[i].c_str() );
        check_file_open( out_files[i], ref_files[i]);

        init_queue[i].set_file_idx( i );

        if(ends_matching)
            ref_start_vector[i].push_back(1);
    }

    std::priority_queue<FileSizeElement> fs_heap(init_queue.begin(), init_queue.end());
    FileSizeElement tmp_fse;
    size_t refn = ref.size();
    for(size_t i = 0; i < refn; ++i)
    {
        tmp_fse = fs_heap.top();
        fs_heap.pop();
        size_t idx = tmp_fse.file_idx;

        tmp_fse.add_file_size( ref[i].seq.length() + 1);
        fs_heap.push( tmp_fse );

        out_files[ idx ] << ">" << i << std::endl;
        for(size_t j = 0; j < ref[i].seq.length(); j += 80)
            out_files[ idx ] << ref[i].seq.substr(j, 80) << std::endl;

        // compute ref_starts
        if(ends_matching)
            ref_start_vector[ idx ].push_back( ref_start_vector[idx].back() + ref[i].seq.length() + 1 );
    }

    // close files
    for(size_t i = 0; i < n; ++i)
        out_files[i].close();
    delete[] out_files;
}

static void save_qry_files(const Multiseq& qry, std::vector<std::string>& qry_files, 
        std::vector<std::string>& rev_qry_files, size_t n)
{
    std::ofstream* out_qry_files, *out_rev_qry_files;
    out_qry_files = new std::ofstream[n];
    out_rev_qry_files = new std::ofstream[n];
    std::vector<FileSizeElement> init_queue(n);

    for(size_t i = 0; i < n; ++i)
    {
        // generate file name
        qry_files[ i ] = get_fname( tmp_path + "/qry.", i, ".fa");
        rev_qry_files[ i ] = get_fname( tmp_path + "/rev_qry.", i, ".fa");
        
        // open files
        out_qry_files[i].open( qry_files[i].c_str() );
        out_rev_qry_files[i].open( rev_qry_files[i].c_str() );
        check_file_open( out_qry_files[i], qry_files[i]);
        check_file_open( out_rev_qry_files[i], rev_qry_files[i]);

        init_queue[ i ].set_file_idx( i );
    }

    std::priority_queue<FileSizeElement> fs_heap(init_queue.begin(), init_queue.end());
    FileSizeElement tmp_fse;
    size_t qryn = qry.size();

    for(size_t qryi = 0; qryi < qryn; ++qryi)
    {
        tmp_fse = fs_heap.top();
        fs_heap.pop();
        size_t idx = tmp_fse.file_idx;

        tmp_fse.add_file_size( qry[ qryi ].seq.length() );
        fs_heap.push( tmp_fse );

        save_fasta(out_qry_files[ idx ], qryi, qry[qryi].seq);
        save_fasta(out_rev_qry_files[ idx ], qryi, qry[qryi].rev);
    }
    
    // close files
    for(size_t i = 0; i < n; ++i)
    {
        out_qry_files[i].close();
        out_rev_qry_files[i].close();
    }
    delete[] out_qry_files;
    delete[] out_rev_qry_files;
}

void mumapi_unsetdir()
{
    tmp_path = "/tmp/mummer_api";
}

std::string mumapi_make_fullpath(const std::string& fname)
{
    return (tmp_path + "/" + fname);
}

void mumapi_set_script_path()
{
    script_path = get_bin_path() + string("/run_nucmer.sh");
}

void mumapi_setdir(const std::string& path)
{
    INFO_PRINT("set tmp path [%s]", path.c_str());
    tmp_path = path;
}

void mumapi_mkdir()
{
    int status;
    status = mkdir(tmp_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if(status == -1 && errno != EEXIST)
    {
        ERROR_PRINT("Error %d: unable to mkdir; %s\n", errno, strerror(errno));
        exit(1);
    }
}

void mumapi_rmdir()
{
#ifndef DONTCLEAN 
    std::string command("rm -rf " + tmp_path);
    system(command.c_str());
#endif
}

void mumapi_setpoolsize(int ps)
{
    if(ps <= 0)
    {
        ERROR_PRINT("poolsize should be at least 1");
        exit(1);
    }
    poolsize = ps;
}

void mumapi_clean_exit(int status)
{
    mumapi_rmdir();
    exit(status);
}

void mumapi_set_engine(int engine)
{
    choose_engine = engine;
}

void make_ref_intervals(const Multiseq& ref, std::vector<size_t>& ref_starts)
{
    if(ref.empty()) return;

    size_t n = ref.size();
    ref_starts.resize( n+1 );
    ref_starts[ 0 ] = 1;
    
    for(size_t i = 0; i < n; ++i)
        ref_starts[i+1] = ref_starts[i] + ref[i].seq.length() + 1;
}

size_t get_interval(const std::vector<size_t>& ref_starts, size_t pos)
{
    if(ref_starts.empty())
        return INF;
    size_t L = 1, R = ref_starts.size() - 2;
    size_t mid;
    while(L <= R)
    {
        mid = (L+R)>>1;
        if(ref_starts[ mid ] <= pos)
            L = mid + 1;
        else
            R = mid - 1;
    }
    return L-1;
}

bool keepthis(const std::vector<size_t>& ref_starts, const std::string& line, size_t qlen)
{
    size_t ref_approx_lpadding, qry_approx_lpadding, matchlen;
    istringstream iss(line);
    iss >> ref_approx_lpadding >> qry_approx_lpadding >> matchlen;

    size_t i = get_interval(ref_starts, ref_approx_lpadding);
    size_t ref_approx_rpadding = ref_starts[i+1] - ref_approx_lpadding - 2;
    size_t qry_approx_rpadding = qlen - qry_approx_lpadding;

    ref_approx_lpadding -= ref_starts[i] - 1;

    return (!((ref_approx_lpadding >= qry_approx_lpadding && 
                ref_approx_rpadding >= qry_approx_rpadding) || 
            (ref_approx_lpadding <= qry_approx_lpadding &&
                ref_approx_rpadding <= qry_approx_rpadding)));
}

void remove_middle_matching(const std::vector<size_t>& ref_starts,
        const Multiseq& qry, 
        const std::string& infile, const std::string& outfile)
{
    std::ifstream fin(infile.c_str());
    check_file_open(fin, infile);
    std::ofstream fout(outfile.c_str());
    check_file_open(fout, outfile);

    std::string line;
    size_t qlen = 0;
    while(getline(fin, line))
    {
        if(line[0] == '>')
        {
            fout << line << endl;
            size_t qryid;
            line[0] = ' ';
            istringstream iss(line);
            iss >> qryid;
            qlen = qry[ qryid ].seq.length();
        }
        else if(keepthis(ref_starts, line, qlen))
            fout << line << endl;
    }

    fin.close();
    fout.close();
}

void exec_external_command_with_pipe(const char* path1, char *const args1[], 
        const char* path2, char *const args2[])
{
    int des_p[2];
    pid_t pid1, pid2;
    if(pipe(des_p) == -1)
    {
        perror("Pipe failed");
        mumapi_clean_exit(1);
    }

    pid1 = fork();
    if(pid1 == 0)
    {
        close(STDOUT_FILENO);
        dup(des_p[1]);
        close(des_p[0]);
        close(des_p[1]);

        execvp(path1, args1);
        perror("execvp error: execvp of the first command failed");
        mumapi_clean_exit(1);
    }
    else if(pid1 < 0)
    {
        perror("The first fork failed");
        mumapi_clean_exit(1);
    }

    pid2 = fork();
    if(pid2 == 0)
    {
        close(STDIN_FILENO);
        dup(des_p[0]);
        close(des_p[1]);
        close(des_p[0]);

        execvp(path2, args2);
        perror("execvp error: execvp of the second command failed");
        mumapi_clean_exit(1);
    }
    else if(pid2 < 0)
    {
        perror("The second fork failed");
        mumapi_clean_exit(1);
    }

    close(des_p[0]);
    close(des_p[1]);

    int status1, status2;
    if(-1 == waitpid(pid1, &status1, 0))
    {
        perror("wait process 1 error");
        mumapi_clean_exit(1);
    }
    if(!WIFEXITED(status1) || WEXITSTATUS(status1) != 0)
    {
        perror("Process 1 failed");
        mumapi_clean_exit(1);
    }

    if(-1 == waitpid(pid2, &status2, 0))
    {
        perror("wait process 2 error");
        mumapi_clean_exit(2);
    }
    if(!WIFEXITED(status2) || WEXITSTATUS(status2) != 0)
    {
        perror("Process 2 failed");
        mumapi_clean_exit(2);
    }
}

void exec_external_command(const char* path, char *const args[])
{
    pid_t pid = fork();
    if(pid == 0)
    {
        execvp(path, args);
        perror("execvp error: execvp failed");
        mumapi_clean_exit(1);
    }
    else if(pid > 1)
    {
        int status;
        if(-1 == waitpid(pid, &status, 0))
        {
            perror("wait process error");
            mumapi_clean_exit(1);
        }
        if(!WIFEXITED(status) || WEXITSTATUS(status) != 0)
        {
            perror("Process failed");
            mumapi_clean_exit(1);
        }
    }
    else
    {
        perror("Fork failed");
        mumapi_clean_exit(1);
    }
}

static void run_prenuc(const MUMTaskParams& params)
{
    /*
    command = std::string("$BIOTOOLS/MUMmer3.23/aux_bin/prenuc ") + ref_file + std::string(" > ") + 
            prefix + std::string(".ntref");
    */
    static const char* bash_path = "/bin/bash";
    char* const args[] = {(char *const)"/bin/bash",
                          (char *const)script_path.c_str(),
                          (char *const)"-P",
                          (char *const)"-r", (char *const)params.ref_file.c_str(),
                          (char *const)"-p", (char *const)params.prefix.c_str(),
                          (char *const)NULL};
    exec_external_command(bash_path, args);
}

static void run_remove_ntref(const MUMTaskParams& params)
{
    static const char* bash_path = "/bin/bash";
    char* const args[] = {(char *const)"bin/bash",
                          (char *const)script_path.c_str(),
                          (char *const)"-F",
                          (char *const)"-p", (char *const)params.prefix.c_str(),
                          (char *const)NULL};
    exec_external_command(bash_path, args);
}

static void run_nucmer(const MUMTaskParams& params, bool forward)
{
    std::string str[6];
    if(forward)
    {
        /*
            oss << "mummer -maxmatch -L -l " << params.min_match_length << " -n " << params.prefix 
                    << ".ntref " << params.qry_file << " 2> " << params.prefix << ".log" <<  " | " 
                << "mgaps -l " << params.min_cluster << " -s " << params.max_gap << " -d " 
                    << params.diagdiff << " -f " << params.diagfactor << " | "
                << "$BIOTOOLS/MUMmer3.23/aux_bin/postnuc -b " << params.breaklen << " " 
                    << params.ref_file << " " << params.qry_file << " " << params.prefix << ".forward";
        */
        static const char* bash_path = "/bin/bash";
        char* const args[] = {(char *const)"/bin/bash",
                              (char *const)script_path.c_str(),
                              (char *const)"-N",
                              (char *const)"-p", (char *const)params.prefix.c_str(),
                              (char *const)"-q", (char *const)params.qry_file.c_str(),
                              (char *const)"-r", (char *const)params.ref_file.c_str(),
                              (char *const)"-l", to_charp( params.min_match_length, str[0] ),
                              (char *const)"-c", to_charp( params.min_cluster, str[1] ),
                              (char *const)"-g", to_charp( params.max_gap, str[2] ),
                              (char *const)"-d", to_charp( params.diagdiff, str[3] ),
                              (char *const)"-f", to_charp( params.diagfactor, str[4]),
                              (char *const)"-b", to_charp( params.breaklen, str[5]),
                              (char *const)NULL};

        exec_external_command(bash_path, args);

    }
    else
    {
        /*
            oss << "mummer -maxmatch -L -l " << params.min_match_length << " -n " << params.prefix 
                    << ".ntref " << params.rev_qry_file << " 2> " << params.prefix << ".log" << " | " 
                << "mgaps -l " << params.min_cluster << " -s " << params.max_gap << " -d " 
                    << params.diagdiff << " -f " << params.diagfactor << " | "
                << "$BIOTOOLS/MUMmer3.23/aux_bin/postnuc -b " << params.breaklen << " " 
                    << params.ref_file << " " << params.rev_qry_file << " " << params.prefix << ".backward"
                << " 2> " << params.prefix << ".log";
        */
        static const char* bash_path = "/bin/bash";
        char* const args[] = {(char *const)"/bin/bash",
                              (char *const)script_path.c_str(),
                              (char *const)"-N",
                              (char *const)"-R",
                              (char *const)"-r", (char *const)params.ref_file.c_str(),
                              (char *const)"-p", (char *const)params.prefix.c_str(),
                              (char *const)"-q", (char *const)params.rev_qry_file.c_str(),
                              (char *const)"-l", to_charp(params.min_match_length, str[0]),
                              (char *const)"-c", to_charp(params.min_cluster, str[1]),
                              (char *const)"-g", to_charp(params.max_gap, str[2]),
                              (char *const)"-d", to_charp(params.diagdiff, str[3]),
                              (char *const)"-f", to_charp(params.diagfactor, str[4]),
                              (char *const)"-b", to_charp(params.breaklen, str[5]),
                              (char *const)NULL};

        exec_external_command(bash_path, args);
    }
}

static void run_mummer(const MUMTaskParams& params, bool forward)
{
    std::string str;
    if(forward)
    {
        /*
        oss << "mummer -maxmatch -L -l " << params.min_match_length << " -n " << params.prefix 
            << ".ntref " << params.qry_file << " > " << params.prefix << ".mums";
        */
        static const char* bash_path = "/bin/bash";
        char* const args[] = {(char *const)"/bin/bash",
                              (char *const)script_path.c_str(),
                              (char *const)"-M",
                              (char *const)"-p", (char *const)params.prefix.c_str(),
                              (char *const)"-q", (char *const)params.qry_file.c_str(),
                              (char *const)"-l", to_charp(params.min_match_length, str),
                              (char *const)NULL};

        exec_external_command(bash_path, args);
    }
    else
    {
        /*
            oss << "mummer -maxmatch -L -l " << params.min_match_length << " -n " << params.prefix 
                    << ".ntref " << params.rev_qry_file << " > " << params.prefix << ".mums";
        */
        static const char* bash_path = "/bin/bash";
        char* const args[] = {(char *const)"/bin/bash",
                              (char *const)script_path.c_str(),
                              (char *const)"-M",
                              (char *const)"-R",
                              (char *const)"-p", (char *const)params.prefix.c_str(),
                              (char *const)"-q", (char *const)params.rev_qry_file.c_str(),
                              (char *const)"-l", to_charp(params.min_match_length, str),
                              (char *const)NULL};

        exec_external_command(bash_path, args);
    }
}

static void run_mgaps_and_postnuc(const MUMTaskParams& params, bool forward)
{
    std::string str[5];
    if(forward)
    {
        /*
            oss << "mgaps -l " << params.min_cluster << " -s " << params.max_gap << " -d " 
                    << params.diagdiff << " -f " << params.diagfactor << " < " << params.prefix << ".lmums | "
                << "$BIOTOOLS/MUMmer3.23/aux_bin/postnuc -b " << params.breaklen << " " << params.ref_file << " " 
                    << params.qry_file << " " << params.prefix << ".forward";
        */
        static const char* bash_path = "/bin/bash";
        char* const args[] = {(char *const)"/bin/bash",
                              (char *const)script_path.c_str(),
                              (char *const)"-C",
                              (char *const)"-r", (char *const)params.ref_file.c_str(),
                              (char *const)"-p", (char *const)params.prefix.c_str(),
                              (char *const)"-q", (char *const)params.qry_file.c_str(),
                              (char *const)"-c", to_charp(params.min_cluster, str[0]),
                              (char *const)"-g", to_charp(params.max_gap, str[1]),
                              (char *const)"-d", to_charp(params.diagdiff, str[2]),
                              (char *const)"-f", to_charp(params.diagfactor, str[3]),
                              (char *const)"-b", to_charp(params.breaklen, str[4]),
                              (char *const)NULL};

        exec_external_command(bash_path, args);
    }
    else
    {
        /*
            oss << "mgaps -l " << params.min_cluster << " -s " << params.max_gap << " -d " 
                    << params.diagdiff << " -f " << params.diagfactor << " < " << params.prefix << ".rmums | "
                << "$BIOTOOLS/MUMmer3.23/aux_bin/postnuc -b " << params.breaklen << " " << params.ref_file << " " 
                    << params.rev_qry_file << " " << params.prefix << ".backward";
        */
        static const char* bash_path = "/bin/bash";
        char* const args[] = {(char *const)"/bin/bash",
                              (char *const)script_path.c_str(),
                              (char *const)"-C",
                              (char *const)"-R",
                              (char *const)"-r", (char *const)params.ref_file.c_str(),
                              (char *const)"-p", (char *const)params.prefix.c_str(),
                              (char *const)"-q", (char *const)params.rev_qry_file.c_str(),
                              (char *const)"-c", to_charp(params.min_cluster, str[0]),
                              (char *const)"-g", to_charp(params.max_gap, str[1]),
                              (char *const)"-d", to_charp(params.diagdiff, str[2]),
                              (char *const)"-f", to_charp(params.diagfactor, str[3]),
                              (char *const)"-b", to_charp(params.breaklen, str[4]),
                              (char *const)NULL};

        exec_external_command(bash_path, args);
    }
}

static void run_delta_filter(const MUMTaskParams& params,bool forward)
{
    if(forward)
    {
        /*
        oss << "delta-filter -g " << params.prefix << ".forward.delta > " << params.prefix << ".forward.filter";
        */
        static const char* bash_path = "/bin/bash";
        char* const args[] = {(char *const)"/bin/bash",
                              (char *const)script_path.c_str(),
                              (char *const)"-D",
                              (char *const)"-p", (char *const)params.prefix.c_str(),
                              (char *const)NULL};
        exec_external_command(bash_path, args);
    }
    else
    {
        /*
        oss << "delta-filter -g " << params.prefix << ".backward.delta > " << params.prefix << ".backward.filter";
        */
        static const char* bash_path = "/bin/bash";
        char* const args[] = {(char *const)"/bin/bash",
                              (char *const)script_path.c_str(),
                              (char *const)"-D",
                              (char *const)"-R",
                              (char *const)"-p", (char *const)params.prefix.c_str(),
                              (char *const)NULL};
        exec_external_command(bash_path, args);
    }
}

static void run_delta_filter(const string& prefix, bool forward)
{
    if(forward)
    {
        static const char* bash_path = "/bin/bash";
        char* const args[] = {(char *const)"/bin/bash",
                              (char *const)script_path.c_str(),
                              (char *const)"-D",
                              (char *const)"-p", (char *const)prefix.c_str(),
                              (char *const)NULL};
        exec_external_command(bash_path, args);
    }
    else
    {
        static const char* bash_path = "/bin/bash";
        char* const args[] = {(char *const)"/bin/bash",
                              (char *const)script_path.c_str(),
                              (char *const)"-D",
                              (char *const)"-R",
                              (char *const)"-p", (char *const)prefix.c_str(),
                              (char *const)NULL};
        exec_external_command(bash_path, args);
    }
}

void run_nucmer_task(void* arg)
{
    //INFO_PRINT("");
    const MUMTaskParams& params = *((MUMTaskParams *)arg);

    if(!params.qry_file.empty() && !params.rev_qry_file.empty())
        run_prenuc(params);

    if(!params.qry_file.empty())
    {
        if(params.remove_middle)
        {
            run_mummer(params, true);
            remove_middle_matching(params.ref_starts, *params.qryp, params.prefix + ".mums", params.prefix + ".lmums");
            run_mgaps_and_postnuc(params, true);
        }
        else
        {
            run_nucmer(params, true);
        }

        run_delta_filter(params, true);
    }

    if(!params.rev_qry_file.empty())
    {
        // backward
        if(params.remove_middle)
        {
            run_mummer(params, false);
            remove_middle_matching(params.ref_starts, *params.qryp, params.prefix + ".mums", params.prefix + ".rmums");
            run_mgaps_and_postnuc(params, false);
        }
        else
        {
            run_nucmer(params, false);
        }

        run_delta_filter(params, false);
    }

    if(!params.qry_file.empty() && !params.rev_qry_file.empty())
        run_remove_ntref(params);
    //INFO_PRINT("");
}

void read_deltas(nucmer_Deltas& deltas, const std::string& fname)
{
    std::ifstream fin(fname.c_str());
    check_file_open(fin, fname);

    std::string line;
    getline(fin, line);
    getline(fin, line);

    while(getline(fin, line))
    {
        if(line[0] == '>')
        {
            deltas.push_back(tigrinc::Nucmer_Delta());
            tigrinc::Nucmer_Delta& t_delta = deltas.back();

            std::istringstream iss(line.substr(1));
            iss >> t_delta.ref_tag >> t_delta.qry_tag >> t_delta.ref_len >> t_delta.qry_len;
            iss.clear();
        }
        else
        {
            deltas.back().clusters.push_back( tigrinc::Nucmer_Cluster() );
            tigrinc::Nucmer_Cluster& t_cluster = deltas.back().clusters.back();

            istringstream iss(line);
            iss >> t_cluster.ref_start >> t_cluster.ref_end >> t_cluster.qry_start >> t_cluster.qry_end >> t_cluster.unknown[0] >> t_cluster.unknown[1] >> t_cluster.unknown[2];
            iss.clear();

            long t;
            while(fin >> t)
            {
                if(t == 0)
                {
                    getline(fin, line);
                    break;
                }
                t_cluster.indels.push_back( t );
            }
        }
    }

    fin.close();
}

void all_v_all_matching(const Multiseq& ref, const Multiseq& qry,
        nucmer_Deltas& deltas_forward, nucmer_Deltas& deltas_backward,
        bool ends_matching /* = false */)
{
    if(choose_engine == MUMMER_ENGINE)
    {
        all_v_all_matching_mummer(ref, qry, deltas_forward, deltas_backward, ends_matching);
    }
    else if(choose_engine == BLASR_ENGINE)
    {
        all_v_all_matching_blasr(ref, qry, deltas_forward, deltas_backward, ends_matching);
    }
    else
    {
        ERROR_PRINT("ERROR: unknown engine: %d", choose_engine);
        mumapi_clean_exit(1);
    }
}

// parallel version
void all_v_all_matching_mummer(const Multiseq& ref, const Multiseq& qry,
        nucmer_Deltas& deltas_forward, nucmer_Deltas& deltas_backward,
        bool ends_matching/* = false */)
{
    static size_t m1 = getmaxtextlenstree();
    static size_t m2 = (Showuint) ~((Uint)0);
    const static size_t t1 = 73;
    size_t n1 = 1, n2;
    unsigned long long S1 = get_total_length(ref) + ref.size();
    unsigned long long S2 = get_total_length(qry);
    long double t2 = (long double)S2 / S1;
    
    size_t k = std::max((size_t)1, (size_t)std::ceil((long double)(S1) * (S2)/((long double)(m1) * m2 * poolsize)));
    if(S1 > 50000000)
        k *= 50;
    n2 = std::max((size_t)1, (size_t)std::ceil( std::max((long double)S2 / m2 * 1.5, std::sqrt(t2 * k * poolsize / t1))));
    if(ref.size() > 1)
        n1 = std::min(std::max((size_t)1, (size_t)std::ceil(k * poolsize / (long double)n2)), ref.size());
    else if(S1 > m1)
    {
        ERROR_PRINT("Reference genome is too large!!!");
        mumapi_clean_exit(1);
    }
    n2 = std::min( std::max((size_t)1, (size_t)std::ceil(std::max((long double)S2 / m2 * 1.5, k * poolsize / (long double)n1))), qry.size() );

    DEBUG_PRINT("S1 = %lld, S2 = %lld, m1 = %lu, m2 = %lu, poolsize = %d, t2 = %Lf", S1, S2, m1, m2, poolsize, t2);
    DEBUG_PRINT("k = %lu, n1 = %lu, n2 = %lu", k, n1, n2);

    std::vector<std::string> ref_files(n1), qry_files(n2), rev_qry_files(n2);
    std::vector<std::vector<size_t> > ref_start_vector(n1);
    
    INFO_PRINT("save ref files");
    save_ref_files(ref, ref_files, ref_start_vector, n1, ends_matching);
    INFO_PRINT("save qry files");
    save_qry_files(qry, qry_files, rev_qry_files, n2);

    INFO_PRINT("create thread pool of size %d, and queue size %lu", poolsize, n1 * n2);
    threadpool_t *pool;
    if((pool = threadpool_create(poolsize, n1 * n2, 0)) == NULL)
    {
        ERROR_PRINT("Create thread pool failed");
        mumapi_clean_exit(1);
    }

    INFO_PRINT("dispatch jobs");
    size_t cnt = 0;
    std::vector<MUMTaskParams* > args;
    args.reserve(n1 * n2);
    for(size_t i1 = 0; i1 < n1; ++i1)
    {
        for(size_t i2 = 0; i2 < n2; ++i2)
        {
            if(ends_matching)
                args.push_back( new MUMTaskParams(ref_files[i1], qry_files[i2],
                                    rev_qry_files[i2], get_fname(tmp_path + "/out.", cnt, ""),
                                    ref_start_vector[i1],
                                    true,
                                    &ref, &qry));
            else
                args.push_back( new MUMTaskParams(ref_files[i1], qry_files[i2],
                                    rev_qry_files[i2], get_fname(tmp_path + "/out.", cnt, "")));
                
            if((threadpool_add(pool, &run_nucmer_task, args[cnt], 0)) != 0)
            {
                ERROR_PRINT("add job (%lu, %lu) failed", i1, i2);
                mumapi_clean_exit(1);
            }
            ++cnt;
        }
    }

    if( (threadpool_destroy(pool, threadpool_graceful)) != 0)
    {
        ERROR_PRINT("Destroy thread pool failed");
        mumapi_clean_exit(1);
    }

    INFO_PRINT("read deltas");

    for(size_t i = 0; i < cnt; ++i)
    {
        read_deltas( deltas_forward, args[i]->prefix + ".forward.filter" );
        read_deltas( deltas_backward, args[i]->prefix + ".backward.filter" );
        delete args[i];
    }
    deltas_forward.trim();
    deltas_backward.trim();
}

/*
classes and functions for using blasr
*/
class Pair
{
public:
    size_t ref_id, qry_id;
public:
    Pair(size_t r_id = 0, size_t q_id = 0):
        ref_id(r_id), qry_id(q_id)
    {
    }

    bool operator<(const Pair& t) const
    {
        return (ref_id < t.ref_id || 
                (ref_id == t.ref_id && qry_id < t.qry_id));
    }
};

class BlasrFormat5Cell
{
public:
    std::string qName;
    long long qLength, qStart, qEnd;
    char qStrand;
    std::string tName;
    long long tLength, tStart, tEnd;
    char tStrand;
    long long score, numMatch, numMismatch, numIns, numDel, mapQV;
    std::string qAlignedSeq, matchPattern, tAlignedSeq;
public:
    bool operator<(const BlasrFormat5Cell& t) const
    {
        return (tStart < t.tStart || 
                (tStart == t.tStart && tEnd < t.tEnd) ||
                (tStart == t.tStart && tEnd == t.tEnd && qStart < t.qStart) ||
                (tStart == t.tStart && tEnd == t.tEnd && qStart == t.qStart && qEnd < t.qEnd));
        
    }
    void process_qName()
    {
        size_t found = qName.rfind('/');
        if(found != std::string::npos)
            qName.erase(found);
    }
};

class BlasrFormat5: public std::vector<BlasrFormat5Cell>
{
public:
    bool operator<(const BlasrFormat5& t) const
    {
        if(this->empty())
        {
            if(t.empty())
                return false;
            else
                return true;
        }
        else
        {
            if(t.empty())
                return false;
            else
            {
                const BlasrFormat5Cell& s_front = this->front();
                const BlasrFormat5Cell& t_front = t.front();
                return (s_front.tName < t_front.tName || 
                        (s_front.tName == t_front.tName && s_front.qName < t_front.qName));
            }
        }
    }

    void save(std::ofstream& fout)
    {
        if(this->empty())
            return;
        const BlasrFormat5Cell& front = this->front();
        fout << '>' << str2lu(front.tName) << ' ' << str2lu( front.qName ) << ' '
            << front.tLength << ' ' << front.qLength << std::endl;
        for(std::vector<BlasrFormat5Cell>::iterator it = this->begin(); it != this->end(); ++it)
        {
            long long score = it->numMismatch + it->numIns + it->numDel;
            fout << (it->tStart + 1) << ' ' << (it->tEnd) << ' '
                << (it->qStart + 1) << ' ' << (it->qEnd) << ' '
                << score << ' ' << score << " 0" << std::endl;

            size_t indels = 1;
            for(size_t i = 0; i < it->matchPattern.length(); ++i)
            {
                if(it->matchPattern[i] == '|')
                    ++indels;
                else
                {
                    if(it->qAlignedSeq[i] == '-')
                    {
                        fout << indels << std::endl;
                        indels = 1;
                    }
                    else if(it->tAlignedSeq[i] == '-')
                    {
                        fout << '-' << indels << std::endl;
                        indels = 1;
                    }
                    else
                        ++indels;
                }
            }
            fout << 0 << std::endl;
        }
    }

    void reverse_complement(std::string& seq)
    {
        for(size_t i = 0; i < seq.length(); ++i)
        {
            switch(seq[i])
            {
                case 'A': case 'a':
                    seq[i] = 'T';
                    break;
                case 'T': case 't':
                    seq[i] = 'A';
                    break;
                case 'G': case 'g':
                    seq[i] = 'C';
                    break;
                case 'C': case 'c':
                    seq[i] = 'G';
                    break;
                case 'N': case 'n': case '-':
                    break;
                default:
                    ERROR_PRINT("Error: unknown nucleotide [%c: %d]", seq[i], int(seq[i]));
                    exit(1);
            }
        }
    }
    
    bool keep_this_cell(BlasrFormat5Cell& cell)
    {
        size_t ref_approx_lpadding = cell.tStart + 1;
        size_t ref_approx_rpadding = cell.tLength - cell.tEnd;
        size_t qry_approx_lpadding = cell.qStart + 1;
        size_t qry_approx_rpadding = cell.qLength - cell.qEnd;

        return (!( (ref_approx_lpadding >= qry_approx_lpadding &&
                        ref_approx_rpadding >= qry_approx_rpadding) || 
                    (ref_approx_lpadding <= qry_approx_lpadding &&
                        ref_approx_rpadding <= qry_approx_rpadding)));
    }

    void push_back(BlasrFormat5Cell& cell, bool ends_matching = false)
    {
        if(cell.tStrand == '-')
        {
            if(cell.qStrand == '+')
            {
                long long tmp_start = cell.qLength - cell.qEnd;
                cell.qEnd = cell.qLength - cell.qStart;
                cell.qStart = tmp_start;
            }
            
            if(ends_matching and not keep_this_cell(cell))
                return;

            std::reverse(cell.qAlignedSeq.begin(), cell.qAlignedSeq.end());
            std::reverse(cell.matchPattern.begin(), cell.matchPattern.end());
            std::reverse(cell.tAlignedSeq.begin(), cell.tAlignedSeq.end());

            reverse_complement(cell.qAlignedSeq);
            reverse_complement(cell.tAlignedSeq);
        }
        else
        {
            if(ends_matching and not keep_this_cell(cell))
                return;
        }
        std::vector<BlasrFormat5Cell>::push_back( cell );
    }
};

static void blasr_save_delta(std::vector<BlasrFormat5>& bf5, std::ofstream& fout)
{
    std::sort(bf5.begin(), bf5.end());
    for(std::vector<BlasrFormat5>::iterator it = bf5.begin(); it != bf5.end(); ++it)
    {
        std::sort(it->begin(), it->end());
        it->save( fout );
    }
}

static void blasr_push_back_cell(std::vector<BlasrFormat5>& bf5, std::map<Pair, size_t>& bf5_indices, BlasrFormat5Cell &cell, bool ends_matching)
{
    Pair pairkey(str2lu(cell.tName), str2lu(cell.qName));
    std::map<Pair, size_t>::iterator finder = bf5_indices.find(pairkey);
    
    size_t idx;
    if(finder == bf5_indices.end())
    {
        idx = bf5.size();
        bf5_indices[ pairkey ] = idx;
        bf5.push_back( BlasrFormat5() );
    }
    else
    {
        idx = finder->second;
    }

    bf5[ idx ].push_back( cell, ends_matching );
}

static void blasr_save_delta_header(std::ofstream& fout, const std::string& refpath, const std::string& qrypath)
{
    fout << refpath << ' ' << qrypath << std::endl;
    fout << "NUCMER" << std::endl;
}

void all_v_all_matching_blasr(const Multiseq& ref, const Multiseq& qry,
        nucmer_Deltas& deltas_forward, nucmer_Deltas& deltas_backward,
        bool ends_matching/* = false */)
{
    std::string ref_file = tmp_path + "/ref.fa";
    std::string qry_file = tmp_path + "/qry.fa";
    std::ofstream fout;

    // save ref
    open_file(fout, ref_file);
    size_t seqn = ref.size();
    for(size_t i = 0; i < seqn; ++i)
        save_fasta(fout, i, ref[i].seq);
    fout.close();

    // save qry
    open_file(fout, qry_file);
    seqn = qry.size();
    for(size_t i = 0; i < seqn; ++i)
        save_fasta(fout, i, qry[i].seq);
    fout.close();

    // run blasr
    blasr_run();

    // transform blasr.out.m5 to blasr.out.forward.delta and blasr.out.backward.delta
    std::ofstream fout_fwd, fout_bwd;
    open_file(fout_fwd, tmp_path + "/blasr.out.forward.delta");
    open_file(fout_bwd, tmp_path + "/blasr.out.backward.delta");
    blasr_save_delta_header(fout_fwd, ref_file, qry_file);
    blasr_save_delta_header(fout_bwd, ref_file, qry_file);

    std::map<Pair, size_t> bf5_fwd_indices, bf5_bwd_indices;
    std::vector<BlasrFormat5> bf5_fwd, bf5_bwd;
    BlasrFormat5Cell cell;
    std::ifstream fin;
    open_file(fin, tmp_path + "/blasr.out.m5");

    while(fin >> cell.qName)
    {
        if(cell.qName == "qName")
            std::getline(fin, cell.qName);
        else
        {
            fin >> cell.qLength >> cell.qStart >> cell.qEnd >> cell.qStrand
                >> cell.tName >> cell.tLength >> cell.tStart >> cell.tEnd >> cell.tStrand
                >> cell.score >> cell.numMatch >> cell.numMismatch >> cell.numIns >> cell.numDel >> cell.mapQV
                >> cell.qAlignedSeq >> cell.matchPattern >> cell.tAlignedSeq;

            cell.process_qName();

            if(cell.qStrand == cell.tStrand)
                blasr_push_back_cell(bf5_fwd, bf5_fwd_indices, cell, ends_matching);
            else
                blasr_push_back_cell(bf5_bwd, bf5_bwd_indices, cell, ends_matching);
            
        }
    }
    fin.close();

    blasr_save_delta(bf5_fwd, fout_fwd);
    blasr_save_delta(bf5_bwd, fout_bwd);

    fout_fwd.close();
    fout_bwd.close();

    // call delta filter on `blasr.out.forward.delta` and `blasr.out.forward.delta`
    run_delta_filter(tmp_path + "/blasr.out", true);
    run_delta_filter(tmp_path + "/blasr.out", false);
    
    // read filter file
    read_deltas(deltas_forward,  tmp_path + "/blasr.out.forward.filter");
    read_deltas(deltas_backward, tmp_path + "/blasr.out.backward.filter");

    deltas_forward.trim();
    deltas_backward.trim();
}

void blasr_align_init(const std::string& X, const std::string& Y)
{
    std::ofstream fout;

    // save ref
    open_file(fout, tmp_path + "/ref.fa");
    save_fasta(fout, 0, X);
    fout.close();

    // save qry
    open_file(fout, tmp_path + "/qry.fa");
    save_fasta(fout, 0, Y);
    fout.close();
}

void blasr_run()
{
    std::string ref_file = tmp_path + "/ref.fa";
    std::string qry_file = tmp_path + "/qry.fa";
    std::string blasr_out = tmp_path + "/blasr.out.m5";
    std::string blasr_unaligned = tmp_path + "/unaligned.fa";
    std::string str[3];
    const static char* blasr_path = "blasr";
    char* const blasr_args[] = {(char *const)"blasr",
                                (char *const)qry_file.c_str(),
                                (char *const)ref_file.c_str(),
                                (char *const)"-m", (char *const)"5",
                                //(char *const)"-out", to_charp(blasr_out, str[0]),
                                (char *const)"-unaligned", to_charp(blasr_unaligned, str[1]),
                                (char *const)"-minMatch", (char *const)"12",
                                (char *const)"-nproc", to_charp(poolsize, str[2]),
                                (char *const)NULL};
    // exec_external_command(blasr_path, blasr_args);
    
    static std::string blasrm5_filter_path = get_bin_path() + "/blasrm5_filter";
    std::string str2[2];
    char* const blasrm5_filter_args[] = {to_charp(blasrm5_filter_path, str2[0]),
                                         (char *const)"-f", (char *const)"20",
                                         (char *const)"-o", to_charp(blasr_out, str2[1]),
                                         (char *const)NULL};
    
    exec_external_command_with_pipe(blasr_path, blasr_args, blasrm5_filter_path.c_str(), blasrm5_filter_args);
}

void blasr_output_transform_forward()
{
    std::ofstream fout_fwd;
    open_file(fout_fwd, tmp_path + "/blasr.out.forward.delta");
    blasr_save_delta_header(fout_fwd, tmp_path + "/ref.fa", tmp_path + "/qry.fa");

    std::map<Pair, size_t> bf5_fwd_indices;
    std::vector<BlasrFormat5> bf5_fwd;
    BlasrFormat5Cell cell;
    std::ifstream fin;
    open_file(fin, tmp_path + "/blasr.out.m5");

    while(fin >> cell.qName)
    {
        if(cell.qName == "qName")
            std::getline(fin, cell.qName);
        else
        {
            fin >> cell.qLength >> cell.qStart >> cell.qEnd >> cell.qStrand
                >> cell.tName >> cell.tLength >> cell.tStart >> cell.tEnd >> cell.tStrand
                >> cell.score >> cell.numMatch >> cell.numMismatch >> cell.numIns >> cell.numDel >> cell.mapQV
                >> cell.qAlignedSeq >> cell.matchPattern >> cell.tAlignedSeq;

            cell.process_qName();

            if(cell.qStrand == cell.tStrand)
                blasr_push_back_cell(bf5_fwd, bf5_fwd_indices, cell, false);
        }
    }
    fin.close();

    blasr_save_delta(bf5_fwd, fout_fwd);
    fout_fwd.close();
}

void blasr_align_choice(const std::string& X, const std::string& Y,
        nucmer_Deltas& deltas, size_t boundary, int choice)
{
    if(choice == 3) // 0b11, align_full
        blasr_align_init(X, Y);
    else if(choice == 1) // 0b01, align_right
        blasr_align_init(X.substr(boundary), Y.length() > X.length() - boundary ? Y.substr(0, X.length() - boundary) : Y);
    else if(choice == 2) // 0b10, align_left
        blasr_align_init(X.substr(0, boundary + 1), Y);
    else
    {
        ERROR_PRINT("Error: unknown choice for doing alignment [%d]", choice);
        mumapi_clean_exit(1);
    }

    blasr_run();
    blasr_output_transform_forward();
    run_delta_filter(tmp_path + "/blasr.out", true);
    read_deltas(deltas, tmp_path + "/blasr.out.forward.filter");
    deltas.trim();
}

// api for mhap
static std::string mhap_script;
class MHAPAln
{// only A_dir == 0 and B_dir == 0 are kept !!!
public:
    size_t A_id, B_id;
    double error, minmers;
    bool A_dir, B_dir;
    size_t A_start, A_end, A_len;
    size_t B_start, B_end, B_len;
public:
    bool keep_alignment_containment(double aligned_ratio, double diff_ratio)
    {
        if(A_dir || B_dir)  return false;
        if(aligned_ratio <= 1e-7 && diff_ratio >= 1 - 1e-7) return true;
        long long A_matched = A_end + 1 - A_start;
        long long B_matched = B_end + 1 - B_start;
        return ((A_matched >= aligned_ratio * A_len || B_matched >= aligned_ratio * B_len) && 
                std::abs(A_matched - B_matched) < diff_ratio * (A_matched + B_matched));
    }

    friend std::istream& operator>>(std::istream& is, MHAPAln& aln)
    {
        if(is >> aln.A_id >> aln.B_id)
        {
            --aln.A_id;     --aln.B_id;
            is >> aln.error >> aln.minmers
                >> aln.A_dir >> aln.A_start >> aln.A_end >> aln.A_len
                >> aln.B_dir >> aln.B_start >> aln.B_end >> aln.B_len;
        }
        return is;
    }
};

class EdgeIdAndWeight
{
public:
    size_t& edge_id;
    long long& weight;
public:
    EdgeIdAndWeight(size_t& e_id, long long& w):
        edge_id(e_id), weight(w)
    {}
};

std::string mhap_get_infile_name(size_t id)
{
    std::ostringstream oss;
    oss << tmp_path << "/mhap_in_" << id << ".fasta";
    return oss.str();
}

std::string mhap_get_outfile_name(size_t id)
{
    std::ostringstream oss;
    oss << tmp_path << "/mhap_out_" << id << ".aln";
    return oss.str();
}

void mhap_pick_best_edge_id(size_t n, std::vector<size_t>& best_edge_id, std::vector<long long>& edge_weight)
{
    mhap_script = get_bin_path() + string("/run_mhap.sh");
    threadpool_t *pool;
    //if((pool = threadpool_create(std::min(poolsize, 8), best_edge_id.size(), 0)) == NULL)
    if((pool = threadpool_create(poolsize, best_edge_id.size(), 0)) == NULL)
    {
        ERROR_PRINT("Create thread pool for mhap failed");
        mumapi_clean_exit(1);
    }

    std::vector<EdgeIdAndWeight *> args(n, NULL);
    for(size_t i = 0; i < n; ++i)
    {
        if(best_edge_id[i] == INF)
        {
            args[ i ] = new EdgeIdAndWeight(best_edge_id[i], edge_weight[i]);
            best_edge_id[ i ] = i; // set this as file ID so as to save space
            if((threadpool_add(pool, &run_mhap_task, args[i], 0)) != 0)
            {
                ERROR_PRINT("add job %lu for mhap failed", i);
                mumapi_clean_exit(1);
            }
        }
    }
    if( (threadpool_destroy(pool, threadpool_graceful)) != 0 )
    {
        ERROR_PRINT("Destroy thread pool for mhap failed");
        mumapi_clean_exit(1);
    }
    for(size_t i = 0; i < n; ++i)
        if(args[i]) delete args[i];
}


void run_mhap_task(void* arg)
{
    EdgeIdAndWeight& ret = *(EdgeIdAndWeight *)arg;
    // ret.edge_id value was set as the file id by the mhap_pick_best_id() function, so as to simplify the thread function call
    std::string infile = mhap_get_infile_name( ret.edge_id );
    std::string outfile = mhap_get_outfile_name( ret.edge_id );

    ret.weight = 0;
    static const char* bash_path = "/bin/bash";
    char* const args[] = {(char *const)"/bin/bash",
                          (char *const)mhap_script.c_str(),
                          (char *const)infile.c_str(),
                          (char *const)outfile.c_str(),
                          (char *const)NULL};
    exec_external_command(bash_path, args);

    // update best_edge_id
    std::ifstream fin;
    open_file(fin, outfile);
    
    size_t res = 0;
    std::vector<size_t> deg;
    std::vector<size_t> deg_bestlen;
    MHAPAln tmp_aln;
    while(fin >> tmp_aln)
    {
        if(tmp_aln.keep_alignment_containment(0.8, 0.2))
        {
            size_t deg_max_size = tmp_aln.A_id > tmp_aln.B_id ? tmp_aln.A_id + 1 : tmp_aln.B_id + 1;
            if(deg_max_size > deg.size())
            {
                deg.resize(deg_max_size, 0);
                deg_bestlen.resize( deg_max_size, 0 );
            }
            // update for A
            ++deg[ tmp_aln.A_id ];
            if(deg_bestlen[ tmp_aln.A_id ] < tmp_aln.A_len)
                deg_bestlen[ tmp_aln.A_id ] = tmp_aln.A_len;
            if(deg[ res ] < deg[ tmp_aln.A_id ] || ( deg[res] == deg[ tmp_aln.A_id ] && deg_bestlen[res] < deg_bestlen[ tmp_aln.A_id ]))
                res = tmp_aln.A_id;
            // update for B
            ++deg[ tmp_aln.B_id ];
            if(deg_bestlen[ tmp_aln.B_id ] < tmp_aln.B_len)
                deg_bestlen[ tmp_aln.B_id ] = tmp_aln.B_len;
            if(deg[ res ] < deg[ tmp_aln.B_id ] || ( deg[res] == deg[ tmp_aln.B_id ] && deg_bestlen[res] < deg_bestlen[ tmp_aln.B_id ]))
                res = tmp_aln.B_id;
        }
    }
    fin.close();

    if(deg.empty())
    {
        ret.edge_id = 0;
        ret.weight = 1;

        size_t best_len = 0;
        open_file(fin, infile);

        char ch;
        size_t fasta_id;
        std::string line;
        while(fin >> ch >> fasta_id)
        {
            fin >> line;
            if(line.length() > best_len)
            {
                best_len = line.length();
                ret.edge_id = fasta_id;
            }
        }
        fin.close();
    }
    else
    {
        ret.edge_id = res;
        ret.weight = deg[ res ];
    }
}

std::string muscle_get_infile_name()
{
    return (tmp_path + "/muscle_in.fasta");
}

std::string muscle_get_outfile_name()
{
    return (tmp_path + "/muscle_out.fasta");
}

void muscle_compute_msa(std::vector<std::vector<size_t> >& nuc_cnt)
{
    static std::string muscle_script = get_bin_path() + string("/run_muscle.sh");
    static std::string infile = muscle_get_infile_name();
    static std::string outfile = muscle_get_outfile_name();

    static const char* bash_path = "/bin/bash";
    char* const args[] = {(char *const)"/bin/bash",
                          (char *const)muscle_script.c_str(),
                          (char *const)infile.c_str(),
                          (char *const)outfile.c_str(),
                          (char *const)NULL};
    exec_external_command(bash_path, args);

    std::ifstream fin;
    open_file(fin, outfile);
    std::string line;

    size_t idx = 0;
    while(std::getline(fin, line))
    {
        if(line[0] == '>')
        {
            idx = 0;
            continue;
        }
        else
        {
            for(size_t i = 0; i < line.length(); ++i, ++idx)
            {
                if(idx >= nuc_cnt.size())
                {
                    if(idx == 0)    nuc_cnt.reserve( line.length() );
                    nuc_cnt.push_back( std::vector<size_t>(5, 0) );
                }
                if(line[i] == 'A')
                    ++nuc_cnt[ idx ][0];
                else if(line[i] == 'T')
                    ++nuc_cnt[ idx ][1];
                else if(line[i] == 'G')
                    ++nuc_cnt[ idx ][2];
                else if(line[i] == 'C')
                    ++nuc_cnt[ idx ][3];
                else
                    ++nuc_cnt[ idx ][4];
            }
        } 
    }
}

} // end of namespace loon
