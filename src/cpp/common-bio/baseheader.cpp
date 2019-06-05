#include <ctime>
#include <cstdarg>
#include "baseheader.h"

#define INFO_PRINT_TEMPLATE \
do{\
    std::time_t raw_time;\
    std::time(&raw_time);\
    struct tm* time_info = std::localtime(&raw_time);\
    fprintf(stderr, "[%02d/%02d/%04d %02d:%02d:%02d]: ",\
            time_info->tm_mon + 1, time_info->tm_mday, time_info->tm_year + 1900,\
            time_info->tm_hour, time_info->tm_min, time_info->tm_sec);\
    va_list arg;\
    va_start(arg, fmt);\
    vfprintf(stderr, fmt, arg);\
    va_end(arg);\
    fprintf(stderr, "\n");\
}while(0)

void INFO_PRINT(const char* fmt, ...)
{
    INFO_PRINT_TEMPLATE;
}

void INFO_PRINT(int N, const char* fmt, ...)
{
    for(int i = 0; i < N; ++i)
        fprintf(stderr, "    ");
    INFO_PRINT_TEMPLATE;
}

