#include "dipole_sources.h"
#include <fstream>

using namespace sloc;
using namespace dealii;

// ----------------------------------------------------------------------------

double DipoleSource::primary_contribution(const Point<3>& field_position)
{
    const double sigma0 = 1.0;
    Point<3> R = field_position - location;
    double R3 = std::pow(R.square(), 1.5);
    return (1.0 / (4 * numbers::PI * sigma0)) * (dipole * R) / R3;
}

// ----------------------------------------------------------------------------

double DipoleSources::primary_contribution(const Point<3>& field_position)
{
    double total = 0;
    for (unsigned int i = 0; i < _sources.size(); ++i)
        total += _sources[i]->primary_contribution(field_position);
    return total;
}

DipoleSources::DipoleSources()
{
}

DipoleSources::~DipoleSources()
{
    clear();
}

void DipoleSources::clear()
{
    std::vector<DipoleSource*>::iterator it;
    for (it = _sources.begin(); it != _sources.end(); ++it)
        delete *it;
    _sources.clear();
}

void DipoleSources::add_source(const Point<3>& location, const Point<3>& dipole)
{
    DipoleSource *src = new DipoleSource;
    src->location = location;
    src->dipole = dipole;
    _sources.push_back(src);
}

void DipoleSources::add_sources(int n, double locations[], double dipoles[])
{
    const double *loc = locations;
    const double *d = dipoles;
    for (int i = 0; i < n; ++i)
    {
        DipoleSource *src = new DipoleSource;
        src->location = Point<3>(loc[3*i+0], loc[3*i+1], loc[3*i+2]);
        src->dipole = Point<3>(d[3*i+0], d[3*i+1], d[3*i+2]);
        _sources.push_back(src);
    }
}

void DipoleSources::read(const char *filename)
{
    unsigned int i;
    unsigned int n_dipoles;
    double x, y, z;
    double px, py, pz;

    clear();

    std::ifstream is;
    is.open(filename);

    is >> n_dipoles;

    if (n_dipoles > 0)
    {
        for (i = 0; i < n_dipoles; ++i)
        {
            // read location
            is >> x;
            is >> y;
            is >> z;

            // read dipole components
            is >> px;
            is >> py;
            is >> pz;

            // make a new source
            DipoleSource *src = new DipoleSource;
            src->location = Point<3>(x, y, z);
            src->dipole = Point<3>(px, py, pz);
            _sources.push_back(src);
        }
    }

    is.close();
}

void DipoleSources::write(const char *filename)
{
    unsigned int i;
    const unsigned int n_dipoles = _sources.size();

    std::ofstream os;
    os.open(filename);

    os << n_dipoles << std::endl;

    if (n_dipoles > 0)
    {
        for (i = 0; i < n_dipoles; ++i)
        {
            DipoleSource *src = _sources[i];
            const Point<3>& location = src->location;
            const Point<3>& dipole = src->dipole;

            os << location(0) << " "
               << location(1) << " "
               << location(2) << " "
               << dipole(0) << " "
               << dipole(1) << " "
               << dipole(2) << " "
               << std::endl;
        }
    }

    os.close();
}

// ----------------------------------------------------------------------------
// EOF
