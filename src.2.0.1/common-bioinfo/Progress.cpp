#include <cmath>
#include <string>
#include "Progress.h"

namespace loon
{

Progress::Progress(): total(1), cur(0), precision(1), accum_precision(0)
{
}

int Progress::back_cnt(double percent = 110)
{
    int ret = float_num + (float_num > 0) + 3;
    if(percent <= 10 - pow(0.1, float_num) + 1e-50)
        return ret + 1;
    else if(percent <= 100 - pow(0.1, float_num) + 1e-50)
        return ret + 2;
    else
        return ret + 3;
}

void Progress::start(double total, int float_num /* = 0*/)
{
    this->total = total;
    this->precision = total * 0.01 * std::pow(0.1, float_num);
    this->cur = 0;
    this->accum_precision = 0;
    this->float_num = float_num;

    std::cerr << "[0%]\b\b\b\b" << std::flush;
}

void Progress::progress(double inc/* = 1*/)
{
    cur += inc;
    accum_precision += inc;
    if(accum_precision >= precision)
    {
        double percent = cur * 100 / total;

        std::streamsize p = std::cerr.precision();

        std::cerr.setf(std::ios::fixed);
        std::cerr.precision( float_num );
        std::cerr << "[" << percent << "%]" << std::string(back_cnt(percent), '\b') << std::flush;
        std::cerr.unsetf(std::ios::fixed);
        std::cerr.precision( p );
        accum_precision = 0;
    }
}

void Progress::stop(bool show100 /* = false */)
{
    if(show100)
        std::cerr << "[100%]" << std::string(float_num + (float_num > 0), ' ') << std::endl;
    else
        std::cerr << std::string(back_cnt(), ' ') << std::string(back_cnt(), '\b') << std::flush;
}

}
