
#include <deal.II/base/point.h>
#include "mesh.h"

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "color_utils.h"

// ----------------------------------------------------------------------------

namespace sloc
{

    dealii::Point<3> centroid(const sloc::Mesh& mesh)
    {
        double c[3] = {0,0,0};
        const int npts = mesh.n_points();
        for (int i = 0; i < npts; i++)
        {
            double x[3];
            mesh.get_point(i,x);

            c[0] += x[0];
            c[1] += x[1];
            c[2] += x[2];
        }
        c[0] /= npts;
        c[1] /= npts;
        c[2] /= npts;
        return dealii::Point<3>(c[0], c[1], c[2]);
    }

}
// ----------------------------------------------------------------------------

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace std;

struct cli_data
{
    fs::path inpath;
    bool debug;
};

void parse_command_line(int *argc, char ***argv, cli_data &d)
{
    try {
        po::options_description desc("Options");
        desc.add_options()
            ("help", "produce help message")
            ("debug", "debug mode")
            ("input,i", po::value<string>(), "input stl file")
            ;
        po::variables_map vm;
        po::store(po::parse_command_line(*argc, *argv, desc), vm);
        po::notify(vm);

        if (vm.count("help") || *argc == 1)
        {
            cout << "Usage: "
                 << ANSI_GREEN << *argv[0] << " [ --debug ] -i INPUTMESH" << ANSI_RESET << endl;
        }

        if (vm.count("input"))
        {
            d.inpath = vm["input"].as<string>();
        }
        else
        {
            cerr << "Error: " << ANSI_RED << "Input argument required!" << ANSI_RESET << endl;
            exit(1);
        }
        
        d.debug = !!vm.count("debug");

    } catch (exception& e) {
        cerr << "Error: " << e.what() << endl;
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

    sloc::Mesh mesh;

    cout << "Reading " << d.inpath.string() << endl;
    mesh.read(d.inpath.c_str());
    cout << "Loaded " << mesh.n_cells() << " cells" << endl;

    dealii::Point<3> center = centroid(mesh);
    cout << "Center: " << center << endl;

    return 0;
}

