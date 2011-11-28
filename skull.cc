#include <iostream>
#include "io_stl.h"

using namespace std;

int main(void)
{
    STL_Mesh mesh;
    stl_read("tmp/Bone2.stl", mesh);
    return 0;
}
