#ifndef IO_GETFEM_H
#define IO_GETFEM_H

#include <string>
#include <vector>
#include <getfem/getfem_mesh.h>

namespace sloc
{

    void read_stl(getfem::mesh& mesh, const char *filename);
    void read_points(std::vector<bgeot::base_node>& points, const char *filename);
    //void write_points(const std::vector<bgeot::base_node>& points, const char *filename);
}


#endif
