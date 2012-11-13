#include <iostream>
#include <deal.II/base/logstream.h>


int main(int argc, const char *argv[])
{
    using namespace std;
    using namespace dealii;

    int level = 3;
    deallog.depth_console(level);
    deallog << argv[0] << endl;

    //cout << "test_dealii" << endl;

    return 0;
}
