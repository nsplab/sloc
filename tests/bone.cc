/* bone.cc
 *
 * Separate Bone2.quad.ucd into two types of interfaces:
 * 
 *   Layer 0 -- interior is bone, exterior is air
 *   Layer 1 -- interior is brain, exterior is bone
 *
 */

#include <deal.II/base/point.h>
#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include "io_ucd.h"

using namespace dealii;
using namespace std;

// ----------------------------------------------------------------------------

struct QuadCell
{
    double verts[4*3];

    void reinit(const std::valarray<double>& nodes, const std::vector<long>& cell_verts)
    {
        // copy first vertex
        verts[3*0 + 0] = nodes[3*cell_verts[0] + 0];
        verts[3*0 + 1] = nodes[3*cell_verts[0] + 1];
        verts[3*0 + 2] = nodes[3*cell_verts[0] + 2];

        // copy second vertex
        verts[3*1 + 0] = nodes[3*cell_verts[1] + 0];
        verts[3*1 + 1] = nodes[3*cell_verts[1] + 1];
        verts[3*1 + 2] = nodes[3*cell_verts[1] + 2];

        // copy third vertex
        verts[3*2 + 0] = nodes[3*cell_verts[2] + 0];
        verts[3*2 + 1] = nodes[3*cell_verts[2] + 1];
        verts[3*2 + 2] = nodes[3*cell_verts[2] + 2];

        // copy fourth vertex
        verts[3*3 + 0] = nodes[3*cell_verts[3] + 0];
        verts[3*3 + 1] = nodes[3*cell_verts[3] + 1];
        verts[3*3 + 2] = nodes[3*cell_verts[3] + 2];
    }
};


void centroid(int npts, double *pts, double c[3])
{
    int i;
    c[0] = c[1] = c[2] = 0.0;
    for (i = 0; i < npts; i++)
    {
        c[0] += pts[3*i+0];
        c[1] += pts[3*i+1];
        c[2] += pts[3*i+2];
    }
    c[0] /= npts;
    c[1] /= npts;
    c[2] /= npts;
}

double distance2(QuadCell& cell)
{
    double c[3];
    centroid(4, cell.verts, c);
    return c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
}

Point<3> cross(Point<3> A, Point<3> B)
{
    //
    //             | i    j   k |
    // C = A x B = | a0  a1  a2 |
    //             | b0  b1  b2 |
    //
    //   = (a1*b2 - a2*b1)i - (a0*b2 - a2*b0)j + (a0*b1 - a1*b0)k
    //
    double C[3];
    C[0] = +(A(1)*B(2) - A(2)*B(1));
    C[1] = -(A(0)*B(2) - A(2)*B(0));
    C[2] = +(A(0)*B(1) - A(1)*B(0));

    return Point<3>(C[0], C[1], C[2]);
}

Point<3> orientation(QuadCell& cell)
{
    //
    // Can determine orientation using triangle ABD
    //
    double *v = cell.verts;
    Point<3> A(v[3*0+0], v[3*0+1], v[3*0+2]);
    Point<3> B(v[3*1+0], v[3*1+1], v[3*1+2]);
    Point<3> D(v[3*3+0], v[3*3+1], v[3*3+2]);
    Point<3> N = cross(B-A, D-A);
    return (N / std::sqrt(N.square()));
}

// ----------------------------------------------------------------------------

int main(void)
{
    string infile = "tmp/Bone2.quad.ucd";
    string outfile = "tmp/Bone2_1.quad.ucd";

    sloc::UCD_File ucd;

    cout << "Reading " << infile << endl;
    ucd.read(infile.c_str());


    //
    // assign material ids
    //
    std::vector<sloc::UCD_Cell*>::iterator it;
    long mat_id;
    QuadCell quad;
    double R2;
    Point<3> cell_normal;
    for (it = ucd._cells.begin(); it != ucd._cells.end(); ++it)
    {
        sloc::UCD_Cell *cell = *it;
        quad.reinit(ucd._nodes, cell->cell_verts);

        R2 = distance2(quad);
        cell_normal = orientation(quad);

        mat_id = 0;

        cell->mat_id = mat_id;
    }

    cout << "Writing " << outfile << endl;
    ucd.write(outfile.c_str());

    return 0;
}

