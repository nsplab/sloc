/* Print known information about mesh file. */

#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc == 1)
    {
        cout << "Usage: " << argv[0] << " <meshfile>" << endl;
        return 0;
    }

    string filename(argv[1]);

    cout << "Opening file '" << filename << "'..."
         << "(functionality not yet implemented)"
         << endl;

    return 0;
}

