
#include "misc_utils.h"

bool is_nan(double x)
{
    // nan is the only value for which the following is true! (although, beware compiler optimizations?)
    // see: http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c
    return (x != x);
}

