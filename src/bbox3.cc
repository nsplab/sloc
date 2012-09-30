#include "bbox3.h"

using namespace std;
using namespace sloc;

std::ostream& operator<<(std::ostream& os, const sloc::bbox3& bbox)
{
    os << "[[" << bbox.min[0] << "," << bbox.max[0] << "]"
       << ",[" << bbox.min[1] << "," << bbox.max[1] << "]"
       << ",[" << bbox.min[2] << "," << bbox.max[2] << "]]";
    return os;
}

void sloc::bbox3::update(double x, double y, double z)
{
    // update x
    if (x < min[0])
        min[0] = x;
    else if (x >= max[0])
        max[0] = x;

    // update y
    if (y < min[1])
        min[1] = y;
    else if (y >= max[1])
        min[1] = y;

    // update z
    if (z < min[2])
        min[2] = z;
    else if (z >= max[2])
        max[2] = z;
}

bool bbox3::within(double x, double y, double z)
{
        return (min[0] <= x) && (x <= max[0]) &&
               (min[1] <= y) && (y <= max[1]) &&
               (min[2] <= z) && (z <= max[2]);
}

bool bbox3::outside(double x, double y, double z)
{
    return !within(x,y,z);
}

