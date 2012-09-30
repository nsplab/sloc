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
#include "io_stl.h"
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

struct cli_data
{
    fs::path inpath;
    fs::path outpath;
    bool debug;
};

void parse_command_line(int *argc, char ***argv, cli_data& d)
{
    try {
        po::options_description desc("Options");
        desc.add_options()
            ("help", "produce help message")
            ("debug", "debug mode")
            ("input,i", po::value<string>(), "input stl file")
            ("output,o", po::value<string>(), "output stl file")
            ;
        po::variables_map vm;
        po::store(po::parse_command_line(*argc, *argv, desc), vm);
        po::notify(vm);

        if (vm.count("help") || *argc == 1)
        {
            cout << "Usage: "
                << ANSI_GREEN << *argv[0] << " [ --debug ] -i inmesh.stl -o outmesh.stl"
                << ANSI_RESET << endl;
            exit(0);
        }

        if (vm.count("input"))
        {
            d.inpath = vm["input"].as<string>();
        }
        else
        {
            cerr << ANSI_RED << "Error: Input argument required!" << ANSI_RESET << endl;
            exit(1);
        }

        if (vm.count("output"))
            d.outpath = vm["output"].as<string>();

        d.debug = !!vm.count("debug");

    } catch (exception& e) {
        cerr << "error: " << e.what() << endl;
        exit(1);
    }
}

int main(int argc, char *argv[])
{
    cli_data d;
    parse_command_line(&argc, &argv, d);

    if (!fs::exists(d.inpath))
    {
        cerr << "File "
             << ANSI_RED << d.inpath.string() << ANSI_RESET
             << " does not exist!" << endl;
        return 1;
    }
    cout << "Reading " << ANSI_GREEN << d.inpath.string() << ANSI_RESET << endl;

    sloc::hptimer *T;

    sloc::STL_File stl_in;
    T = sloc::hp::get_timer("STL_File::read");
    T->start();
    stl_in.read(d.inpath.c_str());
    T->stop();

    const int nel = stl_in.n_facets();
    assert(nel > 0);
    cout << "Loaded " << nel << " facets" << endl;

    valarray<double> area;
    area.resize(nel);

    valarray<double> quality;
    quality.resize(nel);

    T = sloc::hp::get_timer("Facet Loop");
    T->start();
    for (int e = 0; e < nel; e++)
    {
        float x,y,z;

        stl_in.get_facet_vertex(e, 0, x, y, z);
        vcg::Point3<double> a(x,y,z);

        stl_in.get_facet_vertex(e, 1, x, y, z);
        vcg::Point3<double> b(x,y,z);

        stl_in.get_facet_vertex(e, 2, x, y, z);
        vcg::Point3<double> c(x,y,z);

        vcg::Triangle3<double> tri(a,b,c);

        area[e] = vcg::DoubleArea(tri);
        quality[e] = vcg::Quality(a,b,c);
    }
    T->stop();

    if (d.debug)
    {
        savefile("area.dat", area);
        savefile("quality.dat", quality);
    }

    if (!d.outpath.empty())
    {
        // count triangles with non-nan elements
        int nel2 = 0;
        for (int i = 0; i < nel; i++)
        {
            if (!is_nan(area[i]))
                nel2++;
        }

        sloc::STL_File stl_out;
        stl_out.set_facets(nel2);

        int e = 0;
        for (int i = 0; i < nel; i++)
        {
            // skip facets with zero area
            if (is_nan(area[i]))
                continue;

            // copy the facet data
            float nx,ny,nz;
            stl_in.get_normal(i, nx, ny, nz);
            stl_out.set_normal(e, nx, ny, nz);
            for (int j = 0; j < 3; j++)
            {
                float x,y,z;
                stl_in.get_facet_vertex(i, j, x, y, z);
                stl_out.set_facet_vertex(e, j, x, y, z);
            }

            e++;
        }
        cout << "Writing " << ANSI_GREEN << d.outpath.string() << ANSI_RESET << endl;
        stl_out.write(d.outpath.c_str());
        cout << "Wrote " << nel2 << " facets" << endl;
    }

    return 0;
}
