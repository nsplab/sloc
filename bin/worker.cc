// worker.cc - receives jobs and executes them
#include <boost/program_options.hpp>
#include <zmq.hpp>
#include <iostream>

namespace po = boost::program_options;
using namespace std;



int main(int argc, char *argv[])
{
    try
    {
        po::options_description desc("Options");
        desc.add_options()
            ("addr", po::value<string>(), "Address")
            ;
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help"))
        {
            cout << desc << endl;
            return 0;
        }

        cout << "Starting worker..." << endl;

        string addr = "tcp://localhost:5555";
        if (vm.count("addr"))
            addr = vm["addr"].as<string>();

        cout << "Connecting to " << addr << endl;

    }
    catch (std::exception& e)
    {
        cout << e.what() << endl;
    }

    return 0;
}
