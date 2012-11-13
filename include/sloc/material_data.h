#ifndef SLOC_MAT_DATA_H
#define SLOC_MAT_DATA_H

#include <utility>
#include <map>

namespace sloc
{
// ----------------------------------------------------------------------------


class MaterialData
{
public:
    typedef std::pair<double,double> pair_t;
    typedef std::map<unsigned int, pair_t> material_def_t;
    typedef std::map<unsigned int, unsigned int> material_map_t;

public:
    MaterialData();
    ~MaterialData();
    void clear();

public:
    unsigned int n_layers() const;

    void set_layer(unsigned int mat_id, double sigma_int, double sigma_ext);

    double get_sigma_int(unsigned int mat_id) const;
    double get_sigma_ext(unsigned int mat_id) const;
    double get_sigma_avg(unsigned int mat_id) const;

    void read_sigma(const char *filename);
    void write_sigma(const char *filename);

public:
    //
    // XXX: these are temporary APIs (figure out where to place them)
    //

    unsigned int n_cells() const { return _mat.size(); }

    void set_material_id(unsigned int elt_id, unsigned int mat_id);

    int get_material_id(unsigned int elt_id) const;

    void read_materials(const char *filename);
    void write_materials(const char *filename);


public:
    material_def_t _sigma;
    material_map_t _mat;
};

inline unsigned int MaterialData::n_layers() const
{
    return _sigma.size();
}

// ----------------------------------------------------------------------------
}
#endif
