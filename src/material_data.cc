#include "material_data.h"
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
}

double MaterialData::get_sigma_int(unsigned int mat_id) const
{
    const_iterator it = _sigma.find(mat_id);
    if (it != _sigma.end())
    {
        const pair_t& p = (*it).second;
        return p.first;
    }
    return 0;
}

double MaterialData::get_sigma_ext(unsigned int mat_id) const
{
    const_iterator it = _sigma.find(mat_id);
    if (it != _sigma.end())
    {
        const pair_t& p = it->second;
        return p.second;
    }
    return 0;
}

double MaterialData::get_sigma_avg(unsigned int mat_id) const
{
    const_iterator it = _sigma.find(mat_id);
    if (it != _sigma.end())
    {
        const pair_t& p = it->second;
        return (p.first + p.second) / 2;
    }
    return 0;
}

void MaterialData::set_layer(unsigned int mat_id, double sigma_int, double sigma_ext)
{
    _sigma[mat_id] = pair_t(sigma_int, sigma_ext);
}

void MaterialData::read(const char *filename)
{
    unsigned int i;
    unsigned int num_layers;
    unsigned int mat_id;
    double sigma_int, sigma_ext;

    clear();

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

void MaterialData::write(const char *filename)
{

    std::ofstream os;
    os.open(filename);

    const unsigned int num_layers = _sigma.size();
    os << num_layers << endl;

    if (num_layers > 0)
    {
        for (const_iterator it = _sigma.begin(); it != _sigma.end(); ++it)
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
// EOF
