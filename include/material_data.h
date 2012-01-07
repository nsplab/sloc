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
    typedef std::map<unsigned int, pair_t> material_map_t;
    typedef material_map_t::const_iterator const_iterator;

public:
    MaterialData();
    ~MaterialData();

    void clear();

    unsigned int n_layers() const;

    const_iterator begin() const;
    const_iterator end() const;

    double get_sigma_int(unsigned int mat_id) const;
    double get_sigma_ext(unsigned int mat_id) const;
    double get_sigma_avg(unsigned int mat_id) const;

    void set_layer(unsigned int mat_id, double sigma_int, double sigma_ext);

    void read(const char *filename);
    void write(const char *filename);

public:
    material_map_t _sigma;
};

inline unsigned int MaterialData::n_layers() const
{
    return _sigma.size();
}

inline MaterialData::const_iterator MaterialData::begin() const
{
    return _sigma.begin();
}

inline MaterialData::const_iterator MaterialData::end() const
{
    return _sigma.end();
}

// ----------------------------------------------------------------------------
}
#endif
