#include <iostream>
#include <sloc/utils.h>
#include <vcg/space/point3.h>
#include <vcg/space/box3.h>

using namespace std;

template <typename T>
ostream& operator<<(ostream& os, vcg::Point3<T>& p)
{
    os << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")";
    return os;
}

template <typename T>
ostream& operator<<(ostream& os, vcg::Box3<T>& b)
{
    os << "Box3(min=" << b.min << ", max=" << b.max << ")";
    return os;
}

// ----------------------------------------------------------------------------

void test_points()
{
    vcg::Point3d x(0,0,0);
    vcg::Point3d y(4,5,6);

    vcg::Box3d bbox(x,y);
    _PRINT_VALUE(bbox);
    _PRINT_VALUE(bbox.Volume());

    vcg::Point3d z(1,1,1);
    _PRINT_VALUE(z);
    _PRINT_VALUE(bbox.IsIn(z));

    vcg::Point3d z2(7,7,7);
    _PRINT_VALUE(z2);
    _PRINT_VALUE(bbox.IsIn(z2));
}

// ----------------------------------------------------------------------------

int main(int argc, const char *argv[])
{
    test_points();
    return 0;
}

