#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "types.h"
#include "prenuc.h"

Uchar* prenuc(const std::vector<std::string>& s, Uint& len)
{
    char* ret;
    len = 0;

    if(s.empty())
        return NULL;
    for(std::vector<std::string>::const_iterator it = s.begin(); it != s.end(); ++it)
    {
        len += it->length() + 1;
    }

    ret = (char*)malloc(len * sizeof(char));
    if(ret == NULL)
    {
        fprintf(stderr, "Error: Not enough memory\n");
        exit(-1);
    }

    ret[0] = 0;
    for(std::vector<std::string>::const_iterator it = s.begin(); it != s.end(); ++it)
    {
        strcat(ret, it->c_str());
        strcat(ret, "x");
    }
    return (Uchar *)ret;
}
