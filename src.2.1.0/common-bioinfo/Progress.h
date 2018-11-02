#ifndef _PROGRESS_H
#define _PROGRESS_H

#include <iostream>

namespace loon
{
class Progress
{
private:
    double total;
    double precision;
    double cur;
    double accum_precision;
    int float_num;

    int back_cnt(double percent);
public:
    Progress();
    void start(double total, int float_num = 0);
    void progress(double inc = 1);
    void stop(bool show100 = false);
};

}

#endif
