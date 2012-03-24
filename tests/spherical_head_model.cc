/* 
 * Create a 4-layer spherical head model
 *
 */

#include <iostream>
#include <fstream>
#include "geometry_utils.h"
#include "sphere_surface.h"
#include "material_data.h"

using namespace std;
using namespace sloc;

void make_layer(int id, Mesh& mesh, const char *stl_filename, double radius, int levels)
{
    SphereSurface surf;
    Mesh tmesh;

    // make sphere of specified radius, centered at origin.
    surf.set_radius(radius);
    surf.set_refinement(levels);
    surf.write_stl(stl_filename);
    cout << "Wrote " << stl_filename << endl;

    // load the stl file into tmesh
    tmesh.read_stl(stl_filename);

    // convert triangles to quadrilaterals
    tri2quad(tmesh, mesh);

    // tri2quad subdivides triangles using their edge midpoints,
    // so now we have to make sure that all those newly added
    // vertices are also on the surface of the sphere we've defined
    for (int i = 0; i < mesh.n_points(); i++)
    {
        double pt[3];
        mesh.get_point(i, pt);

        double d = sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]);
        if (d < radius)
        {
            pt[0] *= radius/d;
            pt[1] *= radius/d;
            pt[2] *= radius/d;
            mesh.set_point(i, pt);
        }
    }

    // set the material id of our mesh
    for (int e = 0; e < mesh.n_cells(); e++)
        mesh.set_mat(e, id);

    // done making the layer
    return;
}

void append_mesh(const Mesh& source, Mesh& target, int& points_offset, int& cells_offset)
{
    //
    // copy source's data into target's
    //

    int e,n,i;
    const bool verbose = false;

    if (verbose)
    {
        cout << "points_offset = " << points_offset << endl;
        cout << "cells_offset  = " << cells_offset << endl;
    }

    // copy the points
    for (n = 0; n < source.n_points(); n++)
    {
        double pt[3];
        source.get_point(n, pt);
        target.set_point(points_offset + n, pt);
    }

    // copy the cells
    long *cell = new long[source.n_cell_nodes()];
    for (e = 0; e < source.n_cells(); e++)
    {
        source.get_cell(e, cell);

        for (i = 0; i < source.n_cell_nodes(); i++)
            cell[i] += points_offset;

        target.set_cell(cells_offset + e, cell);
    }
    delete [] cell;

    // copy the material ids
    for (e = 0; e < source.n_cells(); e++)
    {
        int mat_id;
        source.get_mat(e, mat_id);
        target.set_mat(cells_offset + e, mat_id);
    }

    // advance the offsets in preparation for the next call
    points_offset += source.n_points();
    cells_offset += source.n_cells();

    // done
    return;
}

int main(void)
{
    // conductivities in units of S/m
    const double sigma_air = 3e-15;
    const double sigma_scalp = 0.275;
    const double sigma_skull = 0.0132;
    const double sigma_csf = 1.79;
    const double sigma_brain = 0.40;
    const double sigma_plasma = 0.667;

    // define the layers
    MaterialData material_data;
    material_data.set_layer(1, sigma_brain, sigma_csf);
    material_data.set_layer(2, sigma_csf, sigma_skull);
    material_data.set_layer(3, sigma_skull, sigma_scalp);
    material_data.set_layer(4, sigma_scalp, sigma_air);
    material_data.write("sphere4.sigma");
    cout << "Wrote sphere4.sigma" << endl;

    // distances in units of m
    const double scalp_thickness = 7e-3;
    const double skull_thickness = 6.5e-3;
    const double csf_thickness = 1.5e-3;
    const double radius_brain = 70e-3;

    // radius of each layer
    double radius1 = radius_brain;
    double radius2 = radius1 + csf_thickness;
    double radius3 = radius2 + skull_thickness;
    double radius4 = radius3 + scalp_thickness;

    // mesh objects
    Mesh head_mesh;
    Mesh layer1, layer2, layer3, layer4;

    // first, build the layers
    int refinement = 0;
    make_layer(1, layer1, "layer1_brain_csf.stl", radius1, refinement);
    make_layer(2, layer2, "layer2_csf_skull.stl", radius2, refinement);
    make_layer(3, layer3, "layer3_skull_scalp.stl", radius3, refinement);
    make_layer(4, layer4, "layer4_scalp_air.stl", radius4, refinement);

    // total number of points and cells over all layers
    int npts = layer1.n_points() + layer2.n_points() + layer3.n_points() + layer4.n_points();
    int ncells = layer1.n_cells() + layer2.n_cells() + layer3.n_cells() + layer4.n_cells();

    // allocate enough memory in head_mesh to do a straight copy
    head_mesh.init_points(npts, 3);
    head_mesh.init_cells(ncells, 4);

    // finally, merge all layers into the head_mesh
    int points_offset = 0;
    int cells_offset = 0;
    append_mesh(layer1, head_mesh, points_offset, cells_offset);
    append_mesh(layer2, head_mesh, points_offset, cells_offset);
    append_mesh(layer3, head_mesh, points_offset, cells_offset);
    append_mesh(layer4, head_mesh, points_offset, cells_offset);

    // all done! let's write out the resulting mesh
    head_mesh.write_ucd("sphere4.surf.ucd");
    cout << "Wrote sphere4.surf.ucd" << endl;

    return 0;
}
