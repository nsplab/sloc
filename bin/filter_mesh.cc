#include <valarray>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "color_utils.h"
#include "misc_utils.h"
#include "mesh.h"
//#include "io_ucd.h"
#include "hptimer.h"

#include <vcg/space/point3.h>
#include <vcg/space/triangle3.h>

template <typename T>
void savefile(std::string filename, std::valarray<T>& xs)
{
    std::cout << "Writing to " << ANSI_MAGENTA << filename << ANSI_RESET << std::endl;
    std::ofstream out(filename.c_str());
    for (int i = 0; i < xs.size(); i++)
        out << i << " " << xs[i] << std::endl;
    out.close();
}

template <typename T>
void savefile(std::string filename, std::vector<T>& xs)
{
    std::cout << "Writing to " << ANSI_MAGENTA << filename << ANSI_RESET << std::endl;
    std::ofstream out(filename.c_str());
    for (int i = 0; i < xs.size(); i++)
        out << xs[i] << std::endl;
    out.close();
}

// ----------------------------------------------------------------------------

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace std;

void parse_command_line(int *argc, char ***argv, fs::path& inpath, fs::path& outpath, bool& debug)
{
	try {
		po::options_description desc("Options");
		desc.add_options()
			("help", "produce help message")
			("debug", "debug mode")
			("input,i", po::value<string>(), "input ucd file")
			("output,o", po::value<string>(), "output ucd file")
			;
		po::variables_map vm;
		po::store(po::parse_command_line(*argc, *argv, desc), vm);
		po::notify(vm);

		if (vm.count("help"))
		{
			cout << "Usage: "
				 << ANSI_GREEN << argv[0] << " [--debug] -i inmesh.ucd -o outmesh.ucd"
				 << ANSI_RESET << endl;
			exit(0);
		}

		if (vm.count("input"))
		{
			inpath = vm["output"].as<string>();
		}
		else
		{
			cerr << ANSI_RED << "Error: Input argument required!" << ANSI_RESET << endl;
			exit(1);
		}

		if (vm.count("output"))
			outpath = vm["output"].as<string>();

		debug = !!vm.count("debug");

	} catch (exception& e) {
		cerr << "error: " << e.what() << endl;
		exit(1);
	}
}

int main(int argc, char *argv[])
{
	fs::path inpath, outpath;
	bool debug;

	parse_command_line(&argc, &argv, inpath, outpath, debug);

	if (!fs::exists(inpath))
	{
		cerr << "File "
			 << ANSI_RED << inpath.filename() << ANSI_RESET
			 << " does not exist!" << endl;
		return 1;
	}
	cout << "Reading " << ANSI_GREEN << inpath.filename() << ANSI_RESET << endl;


	sloc::hptimer *T;

	sloc::UCD_File ucd_in;

	T = sloc::hp::get_timer("UCD_File::read");
	T->start();
	ucd_in.read(inpath.filename().c_str());
	T->stop();

	const int nel = ucd_in.num_cells;
	assert(nel > 0);
	cout << "Loaded " << nel << " cells" << endl;

	valarray<double> area;
	area.resize(nel);

	T = sloc::hp::get_timer("Cell Loop");
	T->start();
	for (int e = 0; e < nel; e++)
	{
	}
	T->stop();

	return 0;
}
