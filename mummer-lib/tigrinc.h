#ifndef  __TIGRINC_HH
#define  __TIGRINC_HH


#include  <cstdio>
#include  <cstdlib>
#include  <cmath>
#include  <cstring>
#include  <cctype>
#include  <climits>
#include  <cfloat>
#include  <ctime>
#include  <cassert>
#include  <cerrno>
#include  <unistd.h>


namespace tigrinc
{
    const int TRUE = 1;
    const int FALSE = 0;
#ifndef  EXIT_FAILURE
    const int EXIT_FAILURE = -1;
#endif
#ifndef  EXIT_SUCCESS
    const int EXIT_SUCCESS = 0;
#endif

    const int INCR_SIZE = 10000;
    const int SMALL_INIT_SIZE = 100;
    const int INIT_SIZE = 10000;
    const int MAX_LINE = 1024;

    FILE *  File_Open  (const char *, const char *);
    void *  Safe_calloc  (size_t, size_t);
    void *  Safe_malloc  (size_t);
    void *  Safe_realloc  (void *, size_t);
    char  Complement  (char);
    bool CompareIUPAC (char, char);
    int  Read_String  (FILE *, char * &, long int &, char [], int);
    void  Reverse_Complement (char S [], long int Lo, long int Hi);
}

#endif
