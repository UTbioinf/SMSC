#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <unistd.h>

using namespace std;

class Args
{
public:
    char* exec;
    double ratio;
    long long filter;
    char* file;
public:
    Args(double r = 2.0, long long f = 0):
        ratio(r), filter(f), file(NULL)
    {}

    void help_msg(bool exit_now = false, int status = 0)
    {
        cerr << "Usage: " << exec << " [options]" << endl;
        cerr << "Options:" << endl;
        cerr << "\t-r <ratio>       Set the ratio for filtering (default 1.0)" << endl;
        cerr << "\t-f <filter>      MapQV smaller than <filter> will be filtered (default 0)" << endl;
        cerr << "\t-o <filename>    Set the filename for output" << endl;
        cerr << "Notice:" << endl;
        cerr << "\tThe input is from `stdin`, and the default output is to `stdout`" << endl;
    
        if(exit_now)
            exit(status);
    }

    double str2double(const string& str)
    {
        double ret;
        istringstream iss(str);
        iss >> ret;
        return ret;
    }

    long long str2longlong(const string& str)
    {
        long long ret;
        istringstream iss(str);
        iss >> ret;
        return ret;
    }

    void check_ratio()
    {
        if(ratio < 0 || ratio > 1)
        {
            cerr << "Error: ratio should be in the range [0, 1]" << endl;
            help_msg(true, 1);
        }
    }

    void parse_parameters(int argc, char* argv[])
    {
        exec = argv[0];
        int opt;
        while((opt = getopt(argc, argv, "r:f:o:h")) != -1)
        {
            switch(opt)
            {
                case 'r':           
                    srand(time(NULL));
                    ratio = str2double(optarg);
                #ifdef DEBUG
                    cerr << "ratio = " << ratio << endl;
                #endif
                    check_ratio();
                    break;
                case 'f':
                    filter = str2longlong(optarg);
                #ifdef DEBUG
                    cerr << "filter = " << filter << endl;
                #endif
                    break;
                case 'o':
                    file = optarg;
                    break;
                case 'h':
                    help_msg(true);
                    break;
                default:
                    cerr << "[" << __FILE__ << ' ' << __LINE__ <<"]: Error: unknown options: " << char(opt) << endl;
                    help_msg(true, 1);
            }
        }
        if(optind != argc)
        {
            cerr << "Error: Too many arguments" << endl;
            help_msg(true, 1);
        }
    }
}args;

void read_n_str(istream& in, int n, vector<string>& strpool)
{
    strpool.resize( n );
    for(int i = 0; i < n; ++i)
    {
        in >> strpool[ i ];
    }
}


int main(int argc, char* argv[])
{
    args.parse_parameters(argc, argv);

    ofstream fout;
    if(args.file != NULL)
    {
        fout.open(args.file);
    }
    
    ostream& out = args.file ? fout : cout;

    string str;
    long long mapQV;
    vector<string> strpool;
    while(cin >> str)
    {
        if(str == "qName")
        {
            out << str;
            getline(cin, str);
            out << str << endl;
        }
        else
        {
            read_n_str(cin, 14, strpool);
            cin >> mapQV;
            if(mapQV < args.filter)
                getline(cin, str);
            else
            {
                if(args.ratio > 1 || rand() / (RAND_MAX+1.0) < args.ratio)
                {
                    out << str;
                    for(size_t i = 0; i < strpool.size(); ++i)
                        out << ' ' << strpool[ i ];
                    out << ' ' << mapQV;
                    getline(cin, str);
                    out << str << endl;
                }
            }
        }
    }

    if(args.file)
        fout.close();
}
