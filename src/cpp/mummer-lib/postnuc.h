#ifndef _POSTNUC_H
#define _POSTNUC_H


#include <vector>
#include <algorithm>
#include "tigrinc.h"
#include "sw_align.h"

namespace tigrinc
{

using namespace std;

const signed char FORWARD_CHAR = 1;
const signed char REVERSE_CHAR = -1;

enum LineType
//-- The type of input line from <stdin>
{
    NO_LINE, HEADER_LINE, MATCH_LINE
};

struct FastaRecord
//-- The essential data of a sequence
{
    char * Id;                 // the fasta ID header tag
    long int len;              // the length of the sequence
    char * seq;                // the sequence data
};


struct Match
//-- An exact match between two sequences A and B
{
    long int sA, sB, len;      // start coordinate in A, in B and the length
};


struct Cluster
//-- An ordered list of matches between two sequences A and B
{
    bool wasFused;             // have the cluster matches been fused yet?
    signed char dirB;          // the query sequence direction
    //      FORWARD_CHAR or REVERSE_CHAR
    vector<Match> matches;     // the ordered set of matches in the cluster
};


struct Synteny
//-- An ordered list of clusters between two sequences A and B
{
    FastaRecord * AfP;         // a pointer to the reference sequence record
    FastaRecord Bf;            // the query sequence record (w/o the sequence)
    vector<Cluster> clusters;  // the ordered set of clusters between A and B
};


struct Alignment
//-- An alignment object between two sequences A and B
{
    signed char dirB;          // the query sequence direction
    long int sA, sB, eA, eB;   // the start in A, B and the end in A, B
    vector<long int> delta;         // the delta values, with NO zero at the end
    long int deltaApos;        // sum of abs(deltas) - #of negative deltas
    //      trust me, it is a very helpful value
    long int Errors, SimErrors, NonAlphas; // errors, similarity errors, nonalphas
};


struct AscendingClusterSort
//-- For sorting clusters in ascending order of their sA coordinate
{
    bool operator() (const Cluster & pA, const Cluster & pB)
    {
        return ( pA.matches.begin( )->sA < pB.matches.begin( )->sA );
    }
};


//------------------------------------------------- Function Declarations ----//
void addNewAlignment(vector<Alignment> & Alignments,
        vector<Cluster>::iterator Cp,
        vector<Match>::iterator Mp);

bool extendBackward(vector<Alignment> & Alignments,
        vector<Alignment>::iterator CurrAp,
        vector<Alignment>::iterator TargetAp, const char * A, const char * B);

bool extendForward(vector<Alignment>::iterator Ap,
        const char * A, long int targetA,
        const char * B, long int targetB, unsigned int m_o);

void extendClusters(vector<Cluster> & Clusters,
        const FastaRecord * Af, const FastaRecord * Bf, FILE * DeltaFile);

void flushAlignments(vector<Alignment> & Alignments,
        const FastaRecord * Af, const FastaRecord * Bf,
        FILE * DeltaFile);

void flushSyntenys(vector<Synteny> & Syntenys, FILE * ClusterFile);

vector<Cluster>::iterator getForwardTargetCluster(vector<Cluster> & Clusters, 
        vector<Cluster>::iterator CurrCp,
        long int & targetA, long int & targetB);

vector<Alignment>::iterator getReverseTargetAlignment(vector<Alignment> & Alignments,
        vector<Alignment>::iterator CurrAp);

bool isShadowedCluster(vector<Cluster>::iterator CurrCp,
        vector<Alignment> & Alignments, vector<Alignment>::iterator Ap);

void parseAbort(const char *);

void parseDelta(vector<Alignment> & Alignments,
        const FastaRecord * Af, const FastaRecord *Bf);

//void processSyntenys(vector<Synteny> & Syntenys,
//        FastaRecord * Af, long int As,
//        FILE * QryFile, FILE * ClusterFile, FILE * DeltaFile);

inline long int revC(long int Coord, long int Len);

void printHelp(const char *);

void printUsage(const char *);

void validateData(vector<Alignment> Alignments,
        vector<Cluster> Clusters,
        const FastaRecord * Af, const FastaRecord * Bf);
}

const int POSTNUC_SET_b = 1;        // set the alignment break (give-up) length to int
const int POSTNUC_SET_B = 2;        // set the diagonal binding for extension to int
const int POSTNUC_SET_d = 4;        // output only match clusters rather than extended alignments; not implemented!!!
const int POSTNUC_SET_e = 8;        // do not extend alignments outward from clusters
const int POSTNUC_SET_s = 0x10;     // don't remove shadowed alignments, useful for aligning
                                    // a sequence to itself to identify repeats
const int POSTNUC_SET_t = 0X20;     // force alignment to ends of sequence if within -b distance

typedef struct
{
    int b;  
    int B;  
    int options;
}postnuc_Params;

typedef struct
{
    void* res;  // store result in this place
    void (*store_header)(void*, const char*, const char*, long int, long int);
    void (*new_cluster)(void*, long int, long int, long int, long int, long int, long int, long int);
    void (*add_delta)(void*, long int);
}postnuc_ProcessResult;

void postnuc_reset_parameters();
void postnuc_set_parameters(postnuc_Params& params);

// this function will set parameters, if provided
void postnuc_init(postnuc_Params* params = NULL);
void postnuc_set_ref(const char* ref_id, long int len, const char* seq);
void postnuc_clear_ref();
void postnuc_set_qry(const char* qry_id, long int len, const char* seq);
void postnuc_clear_qry();
void postnuc_add_Match(long int sA, long int sB, long int len, bool new_cluster = false, bool is_forward = true);
void postnuc_computeDelta(postnuc_ProcessResult& deltaProcess, bool auto_cleanup = false);

#endif
