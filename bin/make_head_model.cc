/*
 * Creates a file head.surf.mesh suitable for bem_forward_solver
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <getfem/getfem_mesh.h>
#include <sloc/io_stl.h>
#include <sloc/material_data.h>

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

// ----------------------------------------------------------------------------

struct parameters
{
    string scalp_stl;
    string outer_skull_stl;
    string inner_skull_stl;
    string brain_stl;
    string arteries_stl;
    string ventricles_stl;
};

void check_file_exists_or_throw(std::string filename)
{
    if (!fs::exists(filename))
    {
        stringstream ss;
        ss << "Invalid file '" << filename << "'";
        throw std::invalid_argument(ss.str());
    }
}

void add_layer_to_mesh(int id, double sigma_ext, double sigma_int, std::string layer_stl, getfem::mesh& mesh, sloc::MaterialData& material_data)
{
    // skip this procedure if file is empty
    if (layer_stl.empty())
        return;

    // but if the filename is non-empty, check that it exists (throw exception if it doesn't)
    check_file_exists_or_throw(layer_stl);

    // load up the stl file
    sloc::STL_File stl;
    cout << "Reading " << layer_stl << endl;
    stl.read(layer_stl.c_str());

    // loop over stl facets
    for (int e = 0; e < stl.n_facets(); ++e)
    {
        bgeot::size_type cv;
        vector<bgeot::size_type> ind(3);
        ind[0] = mesh.add_point(bgeot::base_node(stl.va[3*e+0], stl.va[3*e+1], stl.va[3*e+2]));
        ind[1] = mesh.add_point(bgeot::base_node(stl.vb[3*e+0], stl.vb[3*e+1], stl.vb[3*e+2]));
        ind[2] = mesh.add_point(bgeot::base_node(stl.vc[3*e+0], stl.vc[3*e+1], stl.vc[3*e+2]));
        cv = mesh.add_convex(bgeot::simplex_geotrans(2,1), ind.begin());
        material_data.set_material_id(cv, id);
    }

    material_data.set_layer(id, sigma_int, sigma_ext);
}

void create_head_surf_mesh(parameters& params)
{
    // conductivities in units of S/m
    const double sigma_air = 3e-15;
    const double sigma_scalp = 0.275;
    const double sigma_skull = 0.0132;
    const double sigma_csf = 1.79;
    const double sigma_brain = 0.40;
    const double sigma_plasma = 0.667;

    // mesh stuff
    getfem::mesh surface_mesh;
    sloc::MaterialData material_data;

    // define the layers (NOTE: if stl file is empty string, the layer is skipped)
    add_layer_to_mesh(0, sigma_air, sigma_scalp, params.scalp_stl, surface_mesh, material_data);            // scalp
    add_layer_to_mesh(1, sigma_scalp, sigma_skull, params.outer_skull_stl, surface_mesh, material_data);    // outer skull
    add_layer_to_mesh(2, sigma_skull, sigma_csf, params.inner_skull_stl, surface_mesh, material_data);      // inner skull
    add_layer_to_mesh(3, sigma_csf, sigma_brain, params.brain_stl, surface_mesh, material_data);            // brain
    add_layer_to_mesh(4, sigma_brain, sigma_plasma, params.arteries_stl, surface_mesh, material_data);      // arteries
    add_layer_to_mesh(5, sigma_brain, sigma_csf, params.ventricles_stl, surface_mesh, material_data);       // ventricles

    // write out the conductivities
    material_data.write_sigma("head.surf.sigma");
    cout << "Wrote head.surf.sigma" << endl;

    // write out the materials to a file
    material_data.write_materials("head.surf.mat");
    cout << "Wrote head.surf.mat" << endl;

    // write out the final mesh!
    surface_mesh.write_to_file("head.surf.mesh");
    cout << "Wrote head.surf.mesh" << endl;

    // XXX: export to vtk?
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    try
    {
        string home = fs::path(getenv("HOME")).string();
        string prefix = home + "/meshes";

        parameters p;
        p.scalp_stl = prefix + "/scalp11.stl";
        p.outer_skull_stl = prefix + "/new/outer-skull-20.stl";
        p.inner_skull_stl = prefix + "/new/inner-skull-20.stl";
        p.brain_stl = prefix + "/new/brain2-20.stl";
        p.arteries_stl = prefix + "/artery11.stl";
        p.ventricles_stl = prefix + "/ventricles11.stl";
        create_head_surf_mesh(p);
    }
    catch (std::exception& e)
    {
        cerr << "Error: " << e.what() << endl;
    }
    catch (...)
    {
        cerr << "Exception of unknown type!" << endl;
    }

    return 0;
}
