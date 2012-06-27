#include <octree/octree.h>
#include <cmath>

int main(int argc, char *argv[])
{
    using namespace nomis80;
    Octree<double> o(4096); // Create 4096x4096x4096 octree containing doubles
    o(1,2,3) = M_PI;        // Put pi in (1,2,3)
    o.erase(1,2,3);         // Erase that node
    return 0;
}
