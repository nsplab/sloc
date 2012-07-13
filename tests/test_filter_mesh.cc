#include <valarray>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <boost/filesystem.hpp>

#include "color_utils.h"
#include "io_stl.h"
#include "hptimer.h"

#include <vcg/space/point3.h>
#include <vcg/space/triangle3.h>


template <typename T>
void save(std::string filename, std::valarray<T>& xs)
{
    std::cout << "Writing to " << ANSI_RED << filename << ANSI_RESET << std::endl;
    std::ofstream out(filename.c_str());
    for (int i = 0; i < xs.size(); i++)
        out << i << " " <<  xs[i] << std::endl;
    out.close();
}

template <typename T>
void save(std::string filename, std::vector<T>& xs)
{
    std::cout << "Writing to " << ANSI_RED << filename << ANSI_RESET << std::endl;
    std::ofstream out(filename.c_str());
    for (int i = 0; i < xs.size(); i++)
        out << xs[i] << std::endl;
    out.close();
}


// ----------------------------------------------------------------------------
namespace fs = boost::filesystem;
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage: " << ANSI_RED << argv[0] << " meshfile.stl"
             << ANSI_RESET << endl;
        return 0;
    }

    fs::path p(argv[1]);
    if (!fs::exists(p))
    {
        cerr << "File "
             << ANSI_RED << p.filename()
             << ANSI_RESET << " does not exist!" << endl;
        return 1;
    }
    cout << "Reading " << ANSI_GREEN << p.filename() << ANSI_RESET << endl;

    sloc::hptimer *T;

    sloc::STL_File stl_in;
    T = sloc::hp::get_timer("STL_File::read");
    T->start();
    stl_in.read(p.filename().c_str());
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

    if (true)
    {
        save("area.dat", area);
        save("quality.dat", quality);
    }

    if (true)
    {
        double threshold = 0.0002;

        vector<int> small_facets;
        for (int e = 0; e < nel; e++)
        {
            if (area[e] < threshold)
                small_facets.push_back(e);
        }
        save("small_facets.dat", small_facets);

        const int nsmall = small_facets.size();

        sloc::STL_File stl_out;
        stl_out.set_facets(nel - nsmall);

        int e = 0;
        for (int i = 0; i < nel; i++)
        {
            // skip the small facets
            if (area[i] < threshold)
                continue;

            float nx,ny,nz;
            stl_in.get_normal(i, nx, ny, nz);
            stl_out.set_normal(i, nx, ny, nz);

            for (int j = 0; j < 3; j++)
            {
                float x,y,z;
                stl_in.get_facet_vertex(i, j, x, y, z);
                stl_out.set_facet_vertex(i, j, x, y, z);
            }

            e++;
        }

        cout << "Writing " << ANSI_RED << "foo.stl" << ANSI_RESET << endl;
        stl_out.write("foo.stl");
    }

    //cout << "sloc::hp::report" << endl;
    //sloc::hp::report();
    return 0;
}

