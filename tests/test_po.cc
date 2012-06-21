// https://github.com/boost-lib/program_options/blob/master/example/first.cpp

#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>

namespace po = boost::program_options;
using namespace std;

int main(int argc, char *argv[])
{
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("compression", po::value<double>(), "set compression level")
            ;
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << desc << endl;
            return 0;
        }

        if (vm.count("compression")) {
            cout << "Compression level was set to "
                 << vm["compression"].as<double>()
                 << endl;
        } else {
            cout << "Compression level was not set." << endl;
        }
    } catch (exception& e) {
        cerr << "error: " << e.what() << endl;
    } catch (...) {
        cerr << "Exception of unknown type!" << endl;
    }

    return 0;
}
