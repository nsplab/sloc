// https://github.com/boost-lib/filesystem/blob/master/example/tut3.cpp

#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/filesystem.hpp>

using namespace std;
//using namespace boost::filesystem;
namespace fs = boost::filesystem;

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: " << argv[0] << " path\n";
        return 1;
    }

    fs::path p(argv[1]);   // p reads clearer than argv[1] in the following code

    try
    {
        if (fs::exists(p))    // does p actually exist?
        {
            if (fs::is_regular_file(p))        // is p a regular file?
                cout << p << " size is " << fs::file_size(p) << '\n';
            else if (fs::is_directory(p))      // is p a directory?
            {
                cout << p << " is a directory containing:\n";

                std::copy(fs::directory_iterator(p), fs::directory_iterator(),  // directory_iterator::value_type
                          ostream_iterator<fs::directory_entry>(cout, "\n"));   // is directory_entry, which is
                                                                                // converted to a path by the
                                                                                // path stream inserter
            }
            else
                cout << p << " exists, but is neither a regular file nor a directory\n";
        }
        else
            cout << p << " does not exist\n";
    }
    catch (const fs::filesystem_error& ex)
    {
        cout << ex.what() << '\n';
    }

    return 0;
}
