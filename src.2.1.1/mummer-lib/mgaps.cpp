/* 
*  Modified by zijuexiansheng
*  Programmer:  A. Delcher
*        File:  mgaps.c
*
*  This program reads lists of unique matches between a sequence of strings
*  and a reference string.  For each string in the sequence, it clusters
*  the matches together into groups that may represent longer, inexact
*  matches.
*/

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <climits>
#include <cmath>
#include <cassert>
#include "mgaps.h"

namespace tigrinc
{

const int  DEFAULT_FIXED_SEPARATION = 5;
const long int  DEFAULT_MAX_SEPARATION = 1000;
const long int  DEFAULT_MIN_OUTPUT_SCORE = 200;
const double  DEFAULT_SEPARATION_FACTOR = 0.05;


static inline long int  Max  (long int A, long int B)
    //  Return the larger of  A  and  B .
{
    if  (A < B)
        return  B;
    else
        return  A;
}


static int  Fixed_Separation = DEFAULT_FIXED_SEPARATION;
static long int  Max_Separation = DEFAULT_MAX_SEPARATION;
static long int  Min_Output_Score = DEFAULT_MIN_OUTPUT_SCORE;
static double  Separation_Factor = DEFAULT_SEPARATION_FACTOR;
static int  * UF = NULL;
static int UF_size = 0;
static int  Use_Extents = FALSE;
  // If TRUE use end minus start as length of cluster instead of
  // sum of component lengths


static int  By_Start2(const void * A, const void * B);
static int  By_Cluster(const void * A, const void * B);
static void Filter_Matches(Match_t * A, int & N);
static int  Find(int a);
static int  Process_Cluster(Match_t * A, int N, const char *label, mgaps_add_item add_item, void *clusters);
static void Union(int a, int b);

static int  By_Start2(const void * A, const void * B)
    //  Return how  A  and  B  compare if converted to  Match_t
    //  based on  Start2  value.  If  Start2  values are equal use
    //  Start1  values for comparison.

{
    Match_t  * x, * y;

    x = (Match_t *) A;
    y = (Match_t *) B;

    if  (x -> Start2 < y -> Start2)
        return  -1;
    else if  (x -> Start2 > y -> Start2)
        return  1;
    else if  (x -> Start1 < y -> Start1)
        return  -1;
    else if  (x -> Start1 > y -> Start1)
        return  1;
    else
        return  0;
}



static int  By_Cluster(const void * A, const void * B)
    //  Return how  A  and  B  compare if converted to  Match_t
    //  first based on  cluster_id  value, then by  Start2  value,
    //  then by  Start1  value.

{
    Match_t  * x, * y;

    x = (Match_t *) A;
    y = (Match_t *) B;

    if  (x -> cluster_id < y -> cluster_id)
        return -1;
    else if  (x -> cluster_id > y -> cluster_id)
        return  1;
    else if  (x -> Start2 < y -> Start2)
        return  -1;
    else if  (x -> Start2 > y -> Start2)
        return  1;
    else if  (x -> Start1 < y -> Start1)
        return  -1;
    else if  (x -> Start1 > y -> Start1)
        return  1;
    else
        return  0;
}



static void  Filter_Matches(Match_t * A, int & N)
    //  Remove from  A [0 .. (N - 1)]  any matches that are internal to a repeat,
    //  e.g., if seq1 has 27 As and seq2 has 20 then the first and
    //  last matches will be kept, but the 6 matches in the middle will
    //  be eliminated.  Also combine overlapping matches on the same
    //  diagonal.  Pack all remaining matches into the front of  A  and
    //  reduce the value of  N  if any matches are removed.
    //  Matches in  A  *MUST* be sorted by  Start2  value.

{
    int  i, j;

    for  (i = 0;  i < N;  i ++)
        A [i] . Good = TRUE;

    for  (i = 0;  i < N - 1;  i ++)
    {
        int  i_diag, i_end;

        if  (! A [i] . Good)
            continue;

        i_diag = A [i] . Start2 - A [i] . Start1;
        i_end = A [i] . Start2 + A [i] . Len;

        for  (j = i + 1;  j < N && A [j] . Start2 <= i_end;  j ++)
        {
            int  olap;
            int  j_diag;

            assert (A [i] . Start2 <= A [j] . Start2);

            if  (! A [j] . Good)
                continue;

            j_diag = A [j] . Start2 - A [j] . Start1;
            if  (i_diag == j_diag)
            {
                int  j_extent;

                j_extent = A [j] . Len + A [j] . Start2 - A [i] . Start2;
                if  (j_extent > A [i] . Len)
                {
                    A [i] . Len = j_extent;
                    i_end = A [i] . Start2 + j_extent;
                }
                A [j] . Good = FALSE;
            }
            else if  (A [i] . Start1 == A [j] . Start1)
            {
                olap = A [i] . Start2 + A [i] . Len - A [j] . Start2;
                if  (A [i] . Len < A [j] . Len)
                {
                    if  (olap >=  A [i] . Len / 2)
                    {
                        A [i] . Good = FALSE;
                        break;
                    }
                }
                else if  (A [j] . Len < A [i] . Len)
                {
                    if  (olap >= A [j] . Len / 2)
                    {
                        A [j] . Good = FALSE;
                    }
                }
                else
                {
                    if  (olap >= A [i] . Len / 2)
                    {
                        A [j] . Tentative = TRUE;
                        if  (A [i] . Tentative)
                        {
                            A [i] . Good = FALSE;
                            break;
                        }
                    }
                }
            }
            else if  (A [i] . Start2 == A [j] . Start2)
            {
                olap = A [i] . Start1 + A [i] . Len - A [j] . Start1;
                if  (A [i] . Len < A [j] . Len)
                {
                    if  (olap >=  A [i] . Len / 2)
                    {
                        A [i] . Good = FALSE;
                        break;
                    }
                }
                else if  (A [j] . Len < A [i] . Len)
                {
                    if  (olap >= A [j] . Len / 2)
                    {
                        A [j] . Good = FALSE;
                    }
                }
                else
                {
                    if  (olap >= A [i] . Len / 2)
                    {
                        A [j] . Tentative = TRUE;
                        if  (A [i] . Tentative)
                        {
                            A [i] . Good = FALSE;
                            break;
                        }
                    }
                }
            }
        }
    }

    for  (i = j = 0;  i < N;  i ++)
        if  (A [i] . Good)
        {
            if  (i != j)
                A [j] = A [i];
            j ++;
        }
    N = j;

    for  (i = 0;  i < N;  i ++)
        A [i] . Good = FALSE;

    return;
}



static int  Find(int a)
    //  Return the id of the set containing  a  in  UF .
{
    int  i, j, k;

    if  (UF [a] < 0)
        return  a;

    for  (i = a;  UF [i] > 0;  i = UF [i])
        ;
    for  (j = a;  UF [j] != i;  j = k)
    {
        k = UF [j];
        UF [j] = i;
    }

    return  i;
}





static int Process_Cluster(Match_t * A, int N, const char * label, mgaps_add_item add_item, void *clusters)
    //  Process the cluster of matches in  A [0 .. (N - 1)]  and output them
    //  after a line containing  label .  Return the number of clusters
    //  printed.

{
    long int  adj, total, hi, lo, extent, score;
    int  best, prev;
    int  print_ct = 0;
    int  i, j, k;

    do
    {
        for  (i = 0;  i < N;  i ++)
        {
            A [i] . Simple_Score = A [i] . Len;
            A [i] . Simple_Adj = 0;
            A [i] . Simple_From = -1;
            for  (j = 0;  j < i;  j ++)
            {
                long int  Pen;
                long int  Olap, Olap1, Olap2;

                Olap1 = A [j] . Start1 + A [j] . Len - A [i] . Start1;
                Olap = Max (0, Olap1);
                Olap2 = A [j] . Start2 + A [j] . Len - A [i] . Start2;
                Olap = Max (Olap, Olap2);

                // penalize off diagonal matches
                Pen = Olap + std::abs ( (A [i] . Start2 - A [i] . Start1) -
                        (A [j] . Start2 - A [j] . Start1) );

                if  (A [j] . Simple_Score + A [i] . Len - Pen > A [i] . Simple_Score)
                {
                    A [i] . Simple_From = j;
                    A [i] . Simple_Score = A [j] . Simple_Score + A [i] . Len - Pen;
                    A [i] . Simple_Adj = Olap;
                }
            }
        }

        best = 0;
        for  (i = 1;  i < N;  i ++)
            if  (A [i] . Simple_Score > A [best] . Simple_Score)
                best = i;
        total = 0;
        hi = LONG_MIN;
        lo = LONG_MAX;
        for  (i = best;  i >= 0;  i = A [i] . Simple_From)
        {
            A [i] . Good = TRUE;
            total += A [i] . Len;
            if  (A [i] . Start1 + A [i] . Len > hi)
                hi = A [i] . Start1 + A [i] . Len;
            if  (A [i] . Start1 < lo)
                lo = A [i] . Start1;
        }
        extent = hi - lo;

        if  (Use_Extents)
            score = extent;
        else
            score = total;
        if  (score >= Min_Output_Score)
        {
            print_ct ++;
            prev = -1;
            for  (i = 0;  i < N;  i ++)
                if  (A [i] . Good)
                {
                    if  (prev == -1)
                        add_item(clusters, 1, A[i].Start1,
                                A[i].Start2,
                                A[i].Len,
                                0, 0, 0);
                    else
                    {
                        adj = A [i] . Simple_Adj;

                        add_item(clusters, 0, A[i].Start1 + adj,
                                A[i].Start2 + adj,
                                A[i].Len - adj,
                                -adj,
                                A[i].Start1 + adj - A[prev].Start1 - A[prev].Len,
                                A[i].Start2 + adj - A[prev].Start2 - A[prev].Len);
                    }
                    prev = i;
                }
        }

        for  (i = k = 0;  i < N;  i ++)
            if  (! A [i] . Good)
            {
                if  (i != k)
                {
                    A [k] = A [i];
                }
                k ++;
            }
        N = k;
    }  while  (N > 0);

    return  print_ct;
}







static void Union(int a, int b)
    //  Union the sets whose id's are  a  and  b  in  UF .
{
    assert (UF [a] < 0 && UF [b] < 0);

    if  (UF [a] < UF [b])
    {
        UF [a] += UF [b];
        UF [b] = a;
    }
    else
    {
        UF [b] += UF [a];
        UF [a] = b;
    }

    return;
}

static void Process_Matches(Match_t * A, int N, const char * label, mgaps_add_item add_item, void* clusters)
    //  Process matches  A [1 .. N]  and output them after
    //  a line containing  label .

{
    long int  cluster_size, sep;
    int  print_ct = 0;
    int  a, b, i, j;

    if  (N <= 0)
    {
        return;
    }

    //  Use Union-Find to create connected-components based on
    //  separation and similar diagonals between matches

    if(UF_size <= N)
    {
        if(UF_size == 0)
        {
            UF_size = Max(500, N+1);
            UF = (int*)malloc(UF_size * sizeof(int));

            if(UF == NULL)
            {
                fprintf(stderr, "Error: malloc failed, there is not enough memory\n");
                exit(-1);
            }
        }
        else
        {
            UF_size = Max(UF_size<<1, N+1);
            UF = (int*)realloc(UF, UF_size * sizeof(int));
            if(UF == NULL)
            {
                fprintf(stderr, "Error: realloc failed, there is not enough memory\n");
                exit(-1);
            }
        }
    }
    for  (i = 1;  i <= N;  i ++)
        UF [i] = -1;

    qsort (A + 1, N, sizeof (Match_t), By_Start2);

    Filter_Matches (A + 1, N);

    for  (i = 1;  i < N;  i ++)
    {
        long int  i_end = A [i] . Start2 + A [i] . Len;
        long int  i_diag = A [i] . Start2 - A [i] . Start1;

        for  (j = i + 1;  j <= N;  j ++)
        {
            long int  diag_diff;

            sep = A [j] . Start2 - i_end;
            if  (sep > Max_Separation)
                break;

            diag_diff = std::abs ((A [j] . Start2 - A [j] . Start1) - i_diag);
            if  (diag_diff <= Max (Fixed_Separation, long (Separation_Factor * sep)))
            {
                a = Find (i);
                b = Find (j);
                if  (a != b)
                    Union (a, b);
            }
        }
    }

    //  Set the cluster id of each match

    for  (i = 1;  i <= N;  i ++)
        A [i] . cluster_id = Find (i);

    qsort (A + 1, N, sizeof (Match_t), By_Cluster);

    for  (i = 1;  i <= N;  i += cluster_size)
    {

        for  (j = i + 1;  j <= N && A [i] . cluster_id == A [j] . cluster_id;  j ++)
            ;
        cluster_size = j - i;
        print_ct += Process_Cluster (A + i, cluster_size, label, add_item, clusters);
    }
    return;
}

}

void mgaps_Process_Matches(Match_t * A, int N, const char * label, mgaps_add_item add_item, void* clusters)
{
    tigrinc::Process_Matches(A, N, label, add_item, clusters);
}
 
void mgaps_set_parameters(mgaps_Params* params)
{
    if(params->options & MGAPS_SET_D)
        tigrinc::Fixed_Separation = params->d;
    if(params->options & MGAPS_SET_E)
        tigrinc::Use_Extents = tigrinc::TRUE;
    if(params->options & MGAPS_SET_F)
        tigrinc::Separation_Factor = params->f;
    if(params->options & MGAPS_SET_L)
        tigrinc::Min_Output_Score = params->l;
    if(params->options & MGAPS_SET_S)
        tigrinc::Max_Separation = params->s;
}

void mgaps_ends()
{
    free(tigrinc::UF);
    tigrinc::UF_size = 0;
}
