#ifndef PROGRESS_TIMER_H
#define PROGRESS_TIMER_H

#include <ctime>
#include <iostream>
#include <string>

namespace sloc
{
// ----------------------------------------------------------------------------

class ProgressTimer
{
public:
    ProgressTimer();
    ~ProgressTimer();

    std::string header(const char *firstcol);

    ProgressTimer& start(int total);
    ProgressTimer& update(int current);


public:
    int total;
    int current;
    time_t t_0, t_n;
    double elapsed_mins;
    double mins_per_num;
    double num_per_sec;
    double remaining_mins;
    double total_mins;
    double progress;
};


std::ostream& operator<<(std::ostream& os, const ProgressTimer& T);

// ----------------------------------------------------------------------------
}

#endif
