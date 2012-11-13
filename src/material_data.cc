#include "sloc/material_data.h"
#include <fstream>

using namespace std;
using namespace sloc;

// ----------------------------------------------------------------------------

MaterialData::MaterialData()
{
}

MaterialData::~MaterialData()
{
}

void MaterialData::clear()
{
    _sigma.clear();
    _mat.clear();
}

// ----------------------------------------------------------------------------

void MaterialData::set_layer(unsigned int mat_id, double sigma_int, double sigma_ext)
{
    _sigma[mat_id] = pair_t(sigma_int, sigma_ext);
}

double MaterialData::get_sigma_int(unsigned int mat_id) const
{
    material_def_t::const_iterator it = _sigma.find(mat_id);
    if (it != _sigma.end())
    {
        const pair_t& p = (*it).second;
        return p.first;
    }
    return 0;
}

double MaterialData::get_sigma_ext(unsigned int mat_id) const
{
    material_def_t::const_iterator it = _sigma.find(mat_id);
    if (it != _sigma.end())
    {
        const pair_t& p = it->second;
        return p.second;
    }
    return 0;
}

double MaterialData::get_sigma_avg(unsigned int mat_id) const
{
    material_def_t::const_iterator it = _sigma.find(mat_id);
    if (it != _sigma.end())
    {
        const pair_t& p = it->second;
        return (p.first + p.second) / 2;
    }
    return 0;
}

void MaterialData::read_sigma(const char *filename)
{
    unsigned int i;
    unsigned int num_layers;
    unsigned int mat_id;
    double sigma_int, sigma_ext;

    _sigma.clear();

    std::ifstream is;
    is.open(filename);

    is >> num_layers;

    if (num_layers > 0)
    {
        for (i = 0; i < num_layers; ++i)
        {
            is >> mat_id;
            is >> sigma_int;
            is >> sigma_ext;

            _sigma[mat_id] = pair_t(sigma_int, sigma_ext);
        }
    }

    is.close();
}

void MaterialData::write_sigma(const char *filename)
{
    std::ofstream os;
    os.open(filename);

    const unsigned int num_layers = _sigma.size();
    os << num_layers << endl;

    if (num_layers > 0)
    {
        for (material_def_t::const_iterator it = _sigma.begin(); it != _sigma.end(); ++it)
        {
            unsigned int mat_id = it->first;
            const pair_t& p = it->second;
            os << mat_id << " "
               << p.first << " "
               << p.second << endl;
        }
    }

    os.close();
}

// ----------------------------------------------------------------------------

void MaterialData::set_material_id(unsigned int elt_id, unsigned int mat_id)
{
    _mat[elt_id] = mat_id;
}

int MaterialData::get_material_id(unsigned int elt_id) const
{
    material_map_t::const_iterator it = _mat.find(elt_id);
    if (it != _mat.end())
    {
        int mat_id = it->second;
        return mat_id;
    }
    return -1;
}

void MaterialData::read_materials(const char *filename)
{
    unsigned int i;
    unsigned int num_cells;
    unsigned int elt_id, mat_id;

    _mat.clear();

    std::ifstream is;
    is.open(filename);

    is >> num_cells;

    if (num_cells > 0)
    {
        for (i = 0; i < num_cells; ++i)
        {
            is >> elt_id;
            is >> mat_id;
            _mat[elt_id] = mat_id;
        }
    }

    is.close();
}

void MaterialData::write_materials(const char *filename)
{
    std::ofstream os;
    os.open(filename);

    const unsigned int num_cells = n_cells();

    os << num_cells << endl;

    if (num_cells > 0)
    {
        for (material_map_t::const_iterator it = _mat.begin(); it != _mat.end(); ++it)
        {
            unsigned int elt_id = it->first;
            unsigned int mat_id = it->second;
            os << elt_id << " " << mat_id << endl;
        }
    }

    os.close();
}

// ----------------------------------------------------------------------------
// EOF
