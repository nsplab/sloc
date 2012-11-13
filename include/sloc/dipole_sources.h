#ifndef SLOC_DIPOLES_H
#define SLOC_DIPOLES_H

#include <vector>
#include <deal.II/base/point.h>

namespace sloc
{
// ----------------------------------------------------------------------------

struct DipoleSource
{
    // localization parameter
    dealii::Point<3> location;

    // direction and magnitude of the current source dipole
    dealii::Point<3> dipole;

    double primary_contribution(const dealii::Point<3>& field_position);
};

// ----------------------------------------------------------------------------

class DipoleSources
{
public:
    DipoleSources();
    ~DipoleSources();
    void clear();

    double primary_contribution(const dealii::Point<3>& field_position);

    unsigned int n_sources() const;
    const DipoleSource& operator()(int i) const;

    void add_source(const dealii::Point<3>& location, const dealii::Point<3>& dipole);
    void add_sources(int n, double locations[], double dipoles[]);

    void read(const char *filename);
    void write(const char *filename);

public:
    std::vector<DipoleSource*> _sources;
};

inline unsigned int DipoleSources::n_sources() const
{
    return _sources.size();
}

inline const DipoleSource& DipoleSources::operator()(int i) const
{
    return *(_sources[i]);
}

// ----------------------------------------------------------------------------
}

#endif
