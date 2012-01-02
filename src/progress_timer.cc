#include "progress_timer.h"
#include <iomanip>
#include <sstream>
#include <cassert>
using namespace std;
using namespace sloc;

const int width = 9;

ProgressTimer::ProgressTimer()
{
    t_0 = 0;
    total = 0;
    current = 0;
}

ProgressTimer::~ProgressTimer()
{
}

string ProgressTimer::header(const char *firstcol)
{
    ostringstream hdr;
    hdr << setw(width) << right << firstcol << " "
        << setw(width) << right << "rate" << " "
        << setw(width) << right << "mins" << " "
        << setw(width) << right << "eta" << " "
        << setw(width) << right << "total" << " "
        << setw(width+1) << right << "progress" << endl;
    return hdr.str();
}

ProgressTimer& ProgressTimer::start(int tot)
{
    time(&t_0);
    total = tot;
    return update(0);
}

ProgressTimer& ProgressTimer::update(int cur)
{
    time(&t_n);
    current = cur;
    elapsed_mins = (t_n - t_0) / 60.0;
    mins_per_num = elapsed_mins / current;
    num_per_sec = (1.0/60.0) / mins_per_num;
    remaining_mins = (total - current) * mins_per_num;
    total_mins = total * mins_per_num;
    //progress = 100 * elapsed_mins / total_mins;
    progress = (100.0 * cur) / total;
    return *this;
}

std::ostream&
sloc::operator<<(std::ostream& os, const ProgressTimer& T)
{
    os << setw(width) << right << T.current << " "
       << setw(width) << right << T.num_per_sec << " "
       << setw(width) << right << T.elapsed_mins << " "
       << setw(width) << right << T.remaining_mins << " "
       << setw(width) << right << T.total_mins << " "
       << setw(width) << right << T.progress << "%"
       << "          "
          "          "
          "\r"
       << std::flush;
    return os;
}

