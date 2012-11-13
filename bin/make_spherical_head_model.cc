/* 
 * Create a 4-layer spherical head model
 *
 */

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <sloc/mesh.h>
#include <sloc/material_data.h>
#include <sloc/sphere_surface.h>
#include <sloc/io_getfem.h>
#include <sloc/utils.h>

#include <getfem/getfem_mesh.h>

using namespace std;
using namespace sloc;

// ----------------------------------------------------------------------------

void copy_mesh(const getfem::mesh& src, sloc::Mesh& dst)
{
    // copy points
    int n = 0;
    dst.init_points(src.nb_points(), src.dim());
    for (dal::bv_visitor i(src.points_index()); !i.finished(); ++i, ++n)
    {
        double point[3];
        bgeot::base_node p = src.points()[i];
        point[0] = p[0];
        point[1] = p[1];
        point[2] = p[2];
        dst.set_point(n, point);
    }

    // copy cells
    int e = 0;
    dst.init_cells(src.nb_convex(), src.nb_points_of_convex(0));
    assert(src.nb_points_of_convex(0) == 3);
    for (dal::bv_visitor cv(src.convex_index()); !cv.finished(); ++cv, ++e)
    {
        long cell[3];
        const std::vector<bgeot::size_type>& ind = src.ind_points_of_convex(cv);
        cell[0] = ind[0];
        cell[1] = ind[1];
        cell[2] = ind[2];
        dst.set_cell(e, cell);
    }
}

void copy_mesh(const sloc::Mesh& src, getfem::mesh& dst)
{
    const bool verbose = true;

    if (false)
    {
        _PRINT_VALUE(src.n_cells());
        _PRINT_VALUE(src.n_cell_nodes());
        _PRINT_VALUE(src.n_points());
        _PRINT_VALUE(src.n_dim());
    }

    int e,n;

    assert(src.n_dim() == 3);
    assert(src.n_cell_nodes() == 3);

    double point[src.n_dim()];
    long cell[src.n_cell_nodes()];

    bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(2,1);

    for (e = 0; e < src.n_cells(); ++e)
    {
        src.get_cell(e, cell);

        vector<bgeot::size_type> ind;
        for (n = 0; n < src.n_cell_nodes(); ++n)
        {
            src.get_point(cell[n], point);
            ind.push_back(dst.add_point(bgeot::base_node(point[0], point[1], point[2])));
            //cout << "e = " << e << ", cell[" << n << "] = " << cell[n] << endl;
        }
        //_PRINT_VALUE(ind[0]);
        //_PRINT_VALUE(ind[1]);
        //_PRINT_VALUE(ind[2]);
        dst.add_convex(pgt, ind.begin());
    }
}

void copy_materials(const sloc::Mesh& mesh, sloc::MaterialData& material_data)
{
    int e, mat_id;
    for (e = 0; e < mesh.n_cells(); ++e)
    {
        mesh.get_mat(e, mat_id);
        material_data.set_material_id(e, mat_id);
    }
    material_data.write_materials("sphere4.surf.mat");
}

// ----------------------------------------------------------------------------

void make_layer(int id, Mesh& mesh, const char *stl_filename, double radius, int levels)
{
    SphereSurface surf;

    // make sphere of specified radius, centered at origin.
    surf.set_radius(radius);
    surf.set_refinement(levels);
    surf.write_stl(stl_filename);
    cout << "Wrote " << stl_filename << endl;

    // load the stl file into our sloc::Mesh
    getfem::mesh m;
    sloc::read_stl(m, stl_filename);
    copy_mesh(m, mesh);

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

int main(int argc, char *argv[])
{
    // determine the refinement level
    int refinement = 0;
    if (argc > 1)
        refinement = atoi(argv[1]);
    cout << "Using refinement level " << refinement << endl;

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
    material_data.write_sigma("sphere4.sigma");
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
    make_layer(1, layer1, "layer1_brain_csf.stl", radius1, refinement);
    make_layer(2, layer2, "layer2_csf_skull.stl", radius2, refinement);
    make_layer(3, layer3, "layer3_skull_scalp.stl", radius3, refinement);
    make_layer(4, layer4, "layer4_scalp_air.stl", radius4, refinement);

    // total number of points and cells over all layers
    int npts = layer1.n_points() + layer2.n_points() + layer3.n_points() + layer4.n_points();
    int ncells = layer1.n_cells() + layer2.n_cells() + layer3.n_cells() + layer4.n_cells();

    // allocate enough memory in head_mesh to do a straight copy
    head_mesh.init_points(npts, 3);
    head_mesh.init_cells(ncells, 3);

    // finally, merge all layers into the head_mesh
    int points_offset = 0;
    int cells_offset = 0;
    append_mesh(layer1, head_mesh, points_offset, cells_offset);
    append_mesh(layer2, head_mesh, points_offset, cells_offset);
    append_mesh(layer3, head_mesh, points_offset, cells_offset);
    append_mesh(layer4, head_mesh, points_offset, cells_offset);

    // write out the materials to a file
    copy_materials(head_mesh, material_data);
    material_data.write_materials("sphere4.surf.mat");
    cout << "Wrote sphere4.surf.mat" << endl;

    // all done! let's write out the resulting mesh
    getfem::mesh m;
    copy_mesh(head_mesh, m);
    m.write_to_file("sphere4.surf.mesh");
    cout << "Wrote sphere4.surf.mesh" << endl;

    return 0;
}
