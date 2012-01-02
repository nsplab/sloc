#ifndef IO_STL_H
#define IO_STL_H

#include <string>

namespace sloc
{
// ----------------------------------------------------------------------------

class STL_File
{
public:
    STL_File();
    ~STL_File();
    void read(const char *filename);

    void clear();
    void set_facets(int n);
    void set_normal(int i, float nx, float ny, float nz);
    void set_facet_vertex(int i, int j, float x, float y, float z);

    int n_facets() const;

public:
    int nfacets;
    float *normal;
    float *va;
    float *vb;
    float *vc;
};

inline int STL_File::n_facets() const { return nfacets; }

void stl_read(const std::string& filename, STL_File& mesh);


// ----------------------------------------------------------------------------
}
#endif
