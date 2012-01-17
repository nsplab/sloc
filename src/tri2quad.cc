#include "tri2quad.h"
#include "point_cloud.h"
#include "progress_timer.h"
#include "mesh.h"

using namespace std;
using namespace sloc;

void sloc::centroid(double A[3], double B[3], double C[3], double D[3])
{
    // centroid of three points
    D[0] = (A[0] + B[0] + C[0]) / 3;
    D[1] = (A[1] + B[1] + C[1]) / 3;
    D[2] = (A[2] + B[2] + C[2]) / 3;
}

void sloc::midpoint(double A[3], double B[3], double M[3])
{
    // midpoint coordinates given two edge endpoints
    M[0] = (A[0] + B[0]) / 2;
    M[1] = (A[1] + B[1]) / 2;
    M[2] = (A[2] + B[2]) / 2;
}

void sloc::tri2quad(long t[3], long u[4], long q[3*4])
{
    // This method splits each triangle cell into three quadrilateral cells.
    //
    // We know the triangle connectivity ABC.
    // We can also compute centroid D, and edge midpoints EFG.
    //
    //              C
    // 
    //          G           F
    //                D
    // 
    //     A          E              B
    // 
    // 
    // Connectivity of three new quad elements is just:
    //
    //     AEDG
    //     BFDE
    //     CGDF
    //
    // Note that these orderings preserve the counterclockwise
    // orientation of the original cell.
    //

    long A,B,C,D,E,F,G;

    // vertex ids of triangle being split
    A = t[0]; B = t[1]; C = t[2];

    // ids of the four additional points: centroid and three edge midpoints
    D = u[0]; E = u[1]; F = u[2]; G = u[3];

    // node connectivity of the four new quad elements
    q[4*0+0] = A; q[4*0+1] = E; q[4*0+2] = D; q[4*0+3] = G;
    q[4*1+0] = B; q[4*1+1] = F; q[4*1+2] = D; q[4*1+3] = E;
    q[4*2+0] = C; q[4*2+1] = G; q[4*2+2] = D; q[4*2+3] = F;
}

void sloc::tri2quad(Mesh& tmesh, Mesh& qmesh)
{
    //
    // Split every triangle cell in the first mesh
    // into three quadrilaterals and store it in the
    // second mesh.
    //

    cout << "Calling tri2quad()\n";

    time_t t0, t1;
    ProgressTimer timer;

    double A[3];
    double B[3];
    double C[3];
    double D[3];
    double E[3];
    double F[3];
    double G[3];

    // loop indices
    int e, n;

    // auxiliary array for the new node IDs we'll be inserting
    long *ids;

    //
    // The midpoints of the shared triangle edges give us many
    // duplicate points. We use a point cloud to eliminate duplicate
    // vertex ids.
    //
    PointCloud points;
    points.set_tolerance(1e-8);

    // add points (already presumed unique) from tmesh
    // directly to the point cloud
    // XXX: move this into the PointCloud class
    for (n = 0; n < tmesh.n_points(); n++)
    {
        double *p = new double[3];
        tmesh.get_point(n, p);
        points.points.push_back(p);
    }
    points.npts = tmesh.n_points();

    // four new points per triangle cell
    ids = new long[4 * tmesh.n_cells()];

    cout << "  tmesh points = " << tmesh.n_points() << endl;
    cout << "  tmesh cells  = " << tmesh.n_cells() << endl;

    cout << "  Finding triangle centroids and edge midpoints...\n";
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
        points.add(D[0], D[1], D[2], &ids[4*e+0]);
        points.add(E[0], E[1], E[2], &ids[4*e+1]);
        points.add(F[0], F[1], F[2], &ids[4*e+2]);
        points.add(G[0], G[1], G[2], &ids[4*e+3]);

        // update timer
        if (e % 1000 == 0)
            cout << timer.update(e);
    }
    cout << timer.update(tmesh.n_cells()) << endl;

    time(&t1);
    cout << "  Elapsed time: " << ((t1 - t0) / 60.0) << " mins\n";

    //
    // Now, we initialize the quadrilateral mesh
    //

    // get the final count from the point cloud
    qmesh.init_points(points.n_points(), 3);

    // each triangle is split into three quadrilaterals,
    // so we have three times as many cells
    qmesh.init_cells(3 * tmesh.n_cells(), 4);

    // load point cloud into qmesh
    for (n = 0; n < qmesh.n_points(); n++)
    {
        double pt[3];
        points.get_point(n, pt);
        qmesh.set_point(n, pt);
    }

    cout << "  Splitting triangles into quadrilaterals...\n";
    time(&t0);

    // finally, split each triangle into three quads,
    // and load everything into our quad mesh
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

    time(&t1);
    cout << "  Elapsed time: " << ((t1 - t0) / 60.0) << " mins\n";

    cout << "  qmesh points = " << qmesh.n_points() << endl;
    cout << "  qmesh cells  = " << qmesh.n_cells() << endl;

    delete [] ids;
}

