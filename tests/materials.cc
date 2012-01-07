/* materials.cc
 *
 * For testing the code in material_data.{h,cc}.
 *
 */

#include <iostream>
#include <utility>
#include <string>
#include "material_data.h"

using namespace std;

ostream& operator<<(ostream& os, const sloc::MaterialData& material_data)
{
    sloc::MaterialData::const_iterator it;
    for (it = material_data.begin(); it != material_data.end(); ++it)
    {
        unsigned int mat_id = it->first;
        const std::pair<double,double>& p = it->second;
        os << "Layer " << mat_id << ":"
           << " sigma_int = " << p.first
           << " sigma_ext = " << p.second
           << endl;
    }
    return os;
}

int main(void)
{
    sloc::MaterialData material_data;

    if (true)
    {
        //
        // Read from file
        //

        string filename = "tmp/test.sigma";
        cout << "Reading material data from " << filename << endl;

        material_data.read(filename.c_str());

        cout << material_data << endl;
    }

    if (true)
    {
        material_data.clear();

        //
        // Can also call set_layer()
        //
        material_data.set_layer(0, 1.2, 2.3);
        material_data.set_layer(1, 4.5, 6.7);
        material_data.set_layer(2, 8.9, 0.1);

        cout << "Using the set_layer() method" << endl;
        cout << material_data;

        //
        // Write out to a file
        //
        string filename = "tmp/test2.sigma";
        material_data.write(filename.c_str());
        cout << "Wrote out material data to " << filename << endl;
    }

    return 0;
}
