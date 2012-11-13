#include <sloc/dipole_sources.h>
#include <iostream>

using namespace std;

ostream& operator<<(ostream& os, const sloc::DipoleSources& dipoles)
{
    unsigned int i;
    for (i = 0; i < dipoles.n_sources(); ++i)
    {
        const sloc::DipoleSource& src = dipoles(i);
        os << "Dipole " << i << ": "
           << "Components ("
           << src.dipole(0) << ", "
           << src.dipole(1) << ", "
           << src.dipole(2) << "); "
           << "Location ("
           << src.location(0) << ", "
           << src.location(1) << ", "
           << src.location(2) << ")"
           << endl;
    }
    return os;
}

int main()
{
    sloc::DipoleSources dipoles;

    if (true)
    {
        //
        // Read dipole sources from a file
        //
        string filename = "tmp/test.dipoles";
        dipoles.read(filename.c_str());

        cout << "From file " << filename
            << " -> found " << dipoles.n_sources()
            << " dipoles" << endl;

        cout << dipoles << endl;
    }

    if (true)
    {
        //
        // Set dipole sources from array data
        //
        double dipole_locations[4*3] = {
            0.5, 0, 0,
            0, 0.5, 0,
            0, 0, 0.5,
            0.5, 0.5, 0.5
        };
        double dipole_components[4*3] = {
            0, 1, 0,
            0, 0, 1,
            1, 0, 0,
            0.57735026, 0.57735026, 0.57735026
        };

        dipoles.clear();
        dipoles.add_sources(4, dipole_locations, dipole_components);

        cout << "From array -> found "
            << dipoles.n_sources()
            << " dipoles" << endl;

        cout << dipoles << endl;

        //
        // Write them out to file
        //
        string filename = "tmp/from_array.dipoles";
        dipoles.write(filename.c_str());
        cout << "Wrote out " << filename << endl;
    }

    return 0;
}
