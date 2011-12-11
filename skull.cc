/* Program to convert Bone2.stl and Artery2.stl into .ucd files
 * In the process, we also need to split up every triangle into
 * three quadrilateral elements.
 */
#include <iostream>
#include "io_stl.h"

using namespace std;

int main(void)
{
    STL_Mesh mesh;
    stl_read("tmp/Bone2.stl", mesh);
    return 0;
}
