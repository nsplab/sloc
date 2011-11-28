#ifndef IO_STL_H
#define IO_STL_H

#include <string>

class STL_Mesh
{
public:
    STL_Mesh();
    ~STL_Mesh();
    void clear();
    void set_facets(int n);
    void set_normal(int i, float nx, float ny, float nz);
    void set_facet_vertex(int i, int j, float x, float y, float z);

public:
    int n_facets;
    float *normal;
    float *va;
    float *vb;
    float *vc;
};

void stl_read(const std::string& filename, STL_Mesh& mesh);

#endif
