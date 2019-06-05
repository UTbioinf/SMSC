/*
The parameters for this program are identical to the original mgaps, except:
    * -C is removed
*/
#ifndef _MGAPS_H
#define _MGAPS_H

#include "tigrinc.h"

#define MGAPS_SET_D 1
#define MGAPS_SET_E 2
#define MGAPS_SET_F 4
#define MGAPS_SET_L 8
#define MGAPS_SET_S 0X10

typedef struct
{
    long int  Start1, Start2, Len;
    long int  Simple_Score;
    long int  Simple_From;
    long int  Simple_Adj;
    int  cluster_id : 30;
    unsigned int  Good : 1;
    unsigned int  Tentative : 1;
}Match_t;

typedef struct
{
    int d;      /*Fixed diagonal difference to join matches*/
    double f;   /*Fraction of separation for diagonal difference*/
    long int l; /*Minimum length of cluster match*/
    long int s; /*Maximum separation between matches in cluster*/
    unsigned int options; /*bit options for setting these parameters*/
}mgaps_Params;

/*
The parameters:
1) the variable to save the results
2) 1 if need to create a new cluster, 0 otherwise
3 - 8) the items that should be saved
*/
typedef void (*mgaps_add_item)(void *, int, long int, long int, long int, long int, long int, long int);

/*
The parameters:
1) An array of matchings (starting from 1)
    * `Good` is initialized to FALSE
    * `Tentative` is initialized to FALSE
2) size of the array
3) tag line
4) the callback function that add one item into the clusters
5) the variable that save the item
*/
void mgaps_Process_Matches(Match_t *, int, const char*, mgaps_add_item, void *);
void mgaps_set_parameters(mgaps_Params*);
/*
    Call the following function when you will no longer use mgaps. This function will free some memories.
*/
void mgaps_ends(); 

#endif
