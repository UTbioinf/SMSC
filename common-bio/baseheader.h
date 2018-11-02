#ifndef _BASEHEADER_H
#define _BASEHEADER_H

#include <cstdio>
#include <limits>

#ifndef _IN
    #define _IN
#endif

#ifndef _OUT
    #define _OUT
#endif

#define ZERO 1e-7


namespace loon
{
const size_t INF = std::numeric_limits<std::size_t>::max();
}


#define ERROR_PRINT(...)\
{\
    fprintf(stderr, "[%s %d]: ", __FILE__, __LINE__);\
    fprintf(stderr, ""\
        __VA_ARGS__);\
    fprintf(stderr, "\n");\
}

#define ERROR_PRINT_FILL_BLANK(N, ...)\
{\
    int FILL_BLANK_N = (N);\
    for(int FILL_BLANK_i=0; FILL_BLANK_i<FILL_BLANK_N; ++FILL_BLANK_i)\
        fprintf(stderr, "    ");\
    ERROR_PRINT(__VA_ARGS__);\
}




#define ERROR_BREAKLINE()\
    do{ fprintf(stderr, "\n"); } while(0)

#define INFO_BREAKLINE()\
    do{ fprintf(stderr, "\n"); } while(0)


#ifdef DEBUG
    #define DEBUG_PRINT(...)\
        ERROR_PRINT(__VA_ARGS__);
    #define DEBUG_PRINT_FILL_BLANK(N, ...)\
        ERROR_PRINT_FILL_BLANK(N, __VA_ARGS__);
    #define DEBUG_BREAKLINE()\
        do{ fprintf(stderr, "\n"); } while(0)
#else
    #define DEBUG_PRINT(...)
    #define DEBUG_PRINT_FILL_BLANK(N, ...)
    #define DEBUG_BREAKLINE()
#endif

void INFO_PRINT(const char* fmt, ...);
void INFO_PRINT(int N, const char* fmt, ...);

#endif
