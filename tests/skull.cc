/* skull.cc
 *
 * Program to convert Bone2.stl and Artery2.stl into .ucd files
 * In the process, we also need to split up every triangle into
 * three quadrilateral elements.
 *
 */
#include <iostream>
#include <ctime>
#include "mesh.h"
#include "point_cloud.h"
#include "io_stl.h"
#include "progress_timer.h"

using namespace std;
using namespace sloc;


void centroid(double A[3], double B[3], double C[3], double D[3])
{
    D[0] = (A[0] + B[0] + C[0]) / 3;
    D[1] = (A[1] + B[1] + C[1]) / 3;
    D[2] = (A[2] + B[2] + C[2]) / 3;
}

void midpoint(double A[3], double B[3], double M[3])
{
    M[0] = (A[0] + B[0]) / 2;
    M[1] = (A[1] + B[1]) / 2;
    M[2] = (A[2] + B[2]) / 2;
}

void tri2quad(long t[3], long u[4], long q[3*4])
{
    long A,B,C,D,E,F,G;

    // original vertices
    A = t[0]; B = t[1]; C = t[2];

    // additional points: centroid and edge midpoints
    D = u[0]; E = u[1]; F = u[2]; G = u[3];

    // new quad elements
    q[4*0+0] = A; q[4*0+1] = E; q[4*0+2] = D; q[4*0+3] = G;
    q[4*1+0] = B; q[4*1+1] = F; q[4*1+2] = D; q[4*1+3] = E;
    q[4*2+0] = C; q[4*2+1] = G; q[4*2+2] = D; q[4*2+3] = F;
}

void stl2tmesh2qmesh(STL_File &stl, Mesh &tmesh, Mesh &qmesh)
{
    // Get rid of duplicate points in .stl mesh,
    // and convert every triangle into a quadrilateral.
    // We do everything in one method so we can reuse the
    // point cloud structure built in the first step.
    cout << "Calling stl2tmesh2qmesh()\n";

    time_t t0, t1;
    ProgressTimer timer;

    // loop indices
    int e,n;

    // auxiliary array of node IDs
    long *ids;

    // our triangular mesh tmesh gets its cells from the STL facets.
    int ncells = stl.n_facets();
    int ndim = 3;

    cout << "stl npts = " << (ncells * 3) << endl;
    cout << "stl facets = " << ncells << endl;


    //
    // we can't initialize tmesh just yet.
    //
    // first, we need to determine the list of unique points.
    // we do this by loading them into a point cloud structure
    // that can eliminate duplicates (given a tolerance).
    //
    // when a point is added to the cloud, it is assigned a new id.
    // if the point is already in the cloud, its previous id is returned.
    //
    double epsilon = 1e-8;
    PointCloud points(epsilon);

    ids = new long[ncells * 3];

    cout << "Building point cloud...\n";
    time(&t0);

    cout << timer.header("facets");
    timer.start(ncells);

    for (e = 0; e < ncells; e++)
    {
        points.add(stl.va[3*e+0], stl.va[3*e+1], stl.va[3*e+2], &ids[3*e+0]);
        points.add(stl.vb[3*e+0], stl.vb[3*e+1], stl.vb[3*e+2], &ids[3*e+1]);
        points.add(stl.vc[3*e+0], stl.vc[3*e+1], stl.vc[3*e+2], &ids[3*e+2]);
        if (e % 1000 == 0) cout << timer.update(e);
    }
    cout << timer.update(ncells) << endl;

    time(&t1);
    cout << "Elapsed time: " << ((t1 - t0) / 60.0) << " mins\n";


    // now, we can query the point cloud for its count,
    // and use that to initialize tmesh.
    int npts = points.n_points();

    //
    // initialize our triangular mesh (3 nodes per cell)
    //
    tmesh.init_points(npts, ndim);
    tmesh.init_cells(ncells, 3);

    cout << "tmesh npts = " << tmesh.n_points() << endl;
    cout << "tmesh ncells = " << tmesh.n_cells() << endl;

    // load point cloud into tmesh
    for (n = 0; n < tmesh.n_points(); n++)
    {
        double pt[3];
        points.get_point(n, pt);
        tmesh.set_point(n, pt);
    }

    // load cells into tmesh
    for (e = 0; e < tmesh.n_cells(); e++)
    {
        tmesh.set_cell(e, &ids[3*e]);
    }
    delete [] ids;


    //
    // Now, we split every triangle into three quadrilaterals.
    //

    double A[3];
    double B[3];
    double C[3];
    double D[3];
    double E[3];
    double F[3];
    double G[3];


    // splitting each triangle yields four new points per cell.
    // we reuse the ids array for tracking those extra node IDs.
    ids = new long[4 * tmesh.n_cells()];

    cout << "Finding triangle centroids and edge midpoints...\n";
    time(&t0);

    cout << timer.header("cells");
    timer.start(tmesh.n_cells());

    // calculate new points
    for (e = 0; e < tmesh.n_cells(); e++)
    {
        // load points for triangle e
        long tri_cell[3];
        tmesh.get_cell(e, tri_cell);
        tmesh.get_point(tri_cell[0], A);
        tmesh.get_point(tri_cell[1], B);
        tmesh.get_point(tri_cell[2], C);

        // compute centroid
        centroid(A, B, C, D);

        // midpoint between A and B
        midpoint(A, B, E);

        // midpoint between B and C
        midpoint(B, C, F);

        // midpoint between A and C
        midpoint(A, C, G);

        // add new points to point cloud
        points.add(D[0], D[1], D[3], &ids[4*e+0]);
        points.add(E[0], E[1], E[3], &ids[4*e+1]);
        points.add(F[0], F[1], F[3], &ids[4*e+2]);
        points.add(G[0], G[1], G[3], &ids[4*e+3]);

        // update timer
        if (e % 1000 == 0) cout << timer.update(e);
    }
    cout << timer.update(tmesh.n_cells()) << endl;

    time(&t1);
    cout << "Elapsed time: " << ((t1 - t0) / 60.0) << " mins\n";

    //
    // Now, we initialize our quadrilateral mesh
    //

    // query point cloud to find out how many points we need
    int npts2 = points.n_points();
    qmesh.init_points(npts2, ndim);

    // each triangle is split into three quadrilaterals, so
    // we have three times as many cells
    qmesh.init_cells(3*ncells, 4);

    // load point cloud into qmesh
    for (n = 0; n < npts2; n++)
    {
        double pt[3];
        points.get_point(n, pt);
        qmesh.set_point(n, pt);
    }

    cout << "Splitting triangles into quadrilaterals...\n";
    time(&t0);

    // finally, split each triangle into three quads
    for (e = 0; e < tmesh.n_cells(); e++)
    {
        long tri_cell[3];
        long quad_cells[3*4];
        tmesh.get_cell(e, tri_cell);
        tri2quad(tri_cell, &ids[4*e], quad_cells);
        qmesh.set_cell(3*e+0, &quad_cells[4*0]);
        qmesh.set_cell(3*e+1, &quad_cells[4*1]);
        qmesh.set_cell(3*e+2, &quad_cells[4*2]);
    }
    delete [] ids;

    time(&t1);
    cout << "Elapsed time: " << ((t1 - t0) / 60.0) << " mins\n";

    cout << "qmesh npts = " << qmesh.n_points() << endl;
    cout << "qmesh ncells = " << qmesh.n_cells() << endl;

    // done!
}

int main(void)
{
    STL_File stl;
    Mesh tmesh;
    Mesh qmesh;

    //
    // Process Bone2.stl
    // - read stl file into an object
    // - convert stl into a triangular mesh
    // - convert triangular mesh into quadrilateral mesh
    //
    stl.read("tmp/Bone2.stl");
    stl2tmesh2qmesh(stl, tmesh, qmesh);
    tmesh.write_ucd("tmp/Bone2.tri.ucd");
    qmesh.write_ucd("tmp/Bone2.quad.ucd");

    // do the same for Artery2.stl
    stl.read("tmp/Artery2.stl");
    stl2tmesh2qmesh(stl, tmesh, qmesh);
    tmesh.write_ucd("tmp/Artery2.tri.ucd");
    qmesh.write_ucd("tmp/Bone2.quad.ucd");

    return 0;
}
