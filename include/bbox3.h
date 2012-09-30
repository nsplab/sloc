#ifndef SLOC_BBOX3_H
#define SLOC_BBOX3_H

#include <iostream>

namespace sloc
{
    struct bbox3
    {
        bbox3() {}
        ~bbox3() {}

        void update(double x, double y, double z);
        bool within(double x, double y, double z);
        bool outside(double x, double y, double z);

        double min[3], max[3];
    };

}

std::ostream& operator<<(std::ostream& os, const sloc::bbox3& bbox);

#endif
