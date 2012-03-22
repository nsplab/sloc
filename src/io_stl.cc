#include "io_stl.h"
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <cassert>

using namespace std;
using namespace sloc;

STL_File::STL_File()
{
    nfacets = 0;
    normal = 0;
    va = vb = vc = 0;
}

STL_File::~STL_File()
{
    clear();
}

void STL_File::clear()
{
    if (normal != 0) delete [] normal;
    if (va != 0) delete [] va;
    if (vb != 0) delete [] vb;
    if (vc != 0) delete [] vc;
    nfacets = 0;
    normal = 0;
    va = vb = vb = 0;
}

void STL_File::set_facets(int n)
{
    clear();
    nfacets = n;
    normal = new float[n * 3];
    va = new float[n * 3];
    vb = new float[n * 3];
    vc = new float[n * 3];
}

void STL_File::set_normal(int i, float nx, float ny, float nz)
{
    normal[3*i + 0] = nx;
    normal[3*i + 1] = ny;
    normal[3*i + 2] = nz;
}

void STL_File::set_facet_vertex(int i, int j, float x, float y, float z)
{
    // setting vertex j of triangle i
    if (j == 0) {
        va[3*i + 0] = x;
        va[3*i + 1] = y;
        va[3*i + 2] = z;
    } else if (j == 1) {
        vb[3*i + 0] = x;
        vb[3*i + 1] = y;
        vb[3*i + 2] = z;
    } else if (j == 2) {
        vc[3*i + 0] = x;
        vc[3*i + 1] = y;
        vc[3*i + 2] = z;
    } else {
        // XXX: unreachable! emit error to log...
        assert(false);
    }
}

void STL_File::read(const char *filename)
{
    cout << "Reading " << filename << endl;
    stl_read(filename, *this);
}

void STL_File::write(const char *filename)
{
    cout << "Writing " << filename << endl;
    stl_write(filename, *this);
}

// ----------------------------------------------------------------------------

template <typename T>
T read_binary_type(std::ifstream& in)
{
    char bytes[sizeof(T)];
    in.read(bytes, sizeof(T));
    T* ptr = reinterpret_cast<T*>(bytes);
    return *ptr;
}

void sloc::stl_read(const std::string& filename, STL_File& mesh)
{
    //
    // As described in http://en.wikipedia.org/wiki/STL_(file_format)
    //

    ifstream file;
    unsigned int i,j;
    float x, y, z;
    int attr;

    // open file in binary mode
    file.open(filename.c_str(), ios::in | ios::binary);

    // parse header (80 bytes, ascii)
    char header[80];
    file.read(header, 80);

    // parse number of facets (4 bytes, unsigned long integer)
    unsigned int n = 0;
    n = read_binary_type<uint32_t>(file);

    // initialize the mesh
    mesh.set_facets(n);

    // now, read info for each n facets 
    for (i = 0; i < n; i++)
    {
        // parse normal vector (3 floats - 4 bytes each)
        x = read_binary_type<float>(file);
        y = read_binary_type<float>(file);
        z = read_binary_type<float>(file);
        mesh.set_normal(i, x, y, z);

        // parse vertices (9 floats, 4 bytes each)
        for (j = 0; j < 3; j++)
        {
            x = read_binary_type<float>(file);
            y = read_binary_type<float>(file);
            z = read_binary_type<float>(file);
            mesh.set_facet_vertex(i, j, x, y, z);
        }

        // parse attribute byte count (2 bytes, unsigned integer)
        attr = read_binary_type<uint16_t>(file);
    }

    file.close();
}

void sloc::stl_write(const std::string& filename, STL_File& mesh)
{
    ofstream file;

    // open file in binary mode
    file.open(filename.c_str(), ios::out | ios::binary);

    // first, write the header (zero-filled to use up 80 bytes)
    char header[80];
    memset(&header[0], 0, sizeof(header));
    strcpy(header, "# File created by sloc::STL_File");
    file.write(header, sizeof(header));

    // next, write the number of facets (4 bytes)
    uint32_t n = mesh.n_facets();
    file.write(reinterpret_cast<char*>(&n), sizeof(uint32_t));

    // finally, write out all the facets
    float x, y, z;
    for (int i = 0; i < mesh.n_facets(); i++)
    {
        // write the normal vector (3 floats)
        x = mesh.normal[3*i + 0]; file.write(reinterpret_cast<char*>(&x), sizeof(float));
        y = mesh.normal[3*i + 1]; file.write(reinterpret_cast<char*>(&y), sizeof(float));
        z = mesh.normal[3*i + 2]; file.write(reinterpret_cast<char*>(&z), sizeof(float));

        // write the vertices (3 * 3 floats)
        x = mesh.va[3*i + 0]; file.write(reinterpret_cast<char*>(&x), sizeof(float));
        y = mesh.va[3*i + 1]; file.write(reinterpret_cast<char*>(&y), sizeof(float));
        z = mesh.va[3*i + 2]; file.write(reinterpret_cast<char*>(&z), sizeof(float));

        x = mesh.vb[3*i + 0]; file.write(reinterpret_cast<char*>(&x), sizeof(float));
        y = mesh.vb[3*i + 1]; file.write(reinterpret_cast<char*>(&y), sizeof(float));
        z = mesh.vb[3*i + 2]; file.write(reinterpret_cast<char*>(&z), sizeof(float));

        x = mesh.vc[3*i + 0]; file.write(reinterpret_cast<char*>(&x), sizeof(float));
        y = mesh.vc[3*i + 1]; file.write(reinterpret_cast<char*>(&y), sizeof(float));
        z = mesh.vc[3*i + 2]; file.write(reinterpret_cast<char*>(&z), sizeof(float));

        // write the attribute byte count (2 bytes)
        uint16_t attr = 0;
        file.write(reinterpret_cast<char*>(&attr), sizeof(uint16_t));
    }

    file.close();
}

// EOF
