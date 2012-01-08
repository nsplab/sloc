/* stacked_cubes.cc
 *
 * Create mesh file consisting of two stacked cubes made of two different materials.
 */

#include <vector>
#include <fstream>
#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

using namespace dealii;

void extract_surface(Triangulation<3>& volume_mesh, Triangulation<2,3>& surface_mesh)
{
    if (false)
    {
        // This crashes for some reason... fix this, or use an alternate method
        GridTools::extract_boundary_mesh(volume_mesh, surface_mesh);
    }
    else
    {
        // Try another way...
        // let's build the triangulation directly from a list of cells and vertices

        typedef CellData<2> CellData_2d;
        typedef Point<3> Point_3d;

        std::vector<CellData_2d> cells;
        std::vector<Point_3d> vertices;
        std::vector<bool> touched(volume_mesh.n_vertices(), false);

        // volume vertex indices to surf ones
        std::map<unsigned int, unsigned int> map_vert_index;

        unsigned int v_index;
        CellData_2d c_data;

        // level-0 mesh
        for (Triangulation<3>::cell_iterator cell = volume_mesh.begin(0);
            cell != volume_mesh.end(0);
            ++cell)
        {
            for (unsigned int i = 0; i < (GeometryInfo<2>::faces_per_cell); ++i)
            {
                Triangulation<3>::face_iterator face = cell->face(i);

                if (face->at_boundary())
                {
                    for (unsigned int j = 0; j < (GeometryInfo<2>::vertices_per_cell); ++j)
                    {
                        v_index = face->vertex_index(j);

                        if (!touched[v_index])
                        {
                            vertices.push_back(face->vertex(j));
                            map_vert_index[v_index] = vertices.size() - 1;
                            touched[v_index] = true;
                        }

                        c_data.vertices[j] = map_vert_index[v_index];
                        c_data.material_id = face->boundary_indicator();
                    }
                }

                cells.push_back(c_data);
            }
        }

        Assert(cells.size() > 0, ExcMessage("No boundary faces selected"));

        surface_mesh.create_triangulation(vertices, cells, SubCellData());
    }
}

int main(void)
{
    std::string outfile = "tmp/stacked_cubes.ucd";

    Triangulation<3> volume_mesh;
    Triangulation<2,3> surface_mesh;

    std::vector<unsigned int> n_subdivisions;
    n_subdivisions.push_back(5);
    n_subdivisions.push_back(5);
    n_subdivisions.push_back(10);

    Point<3> bottom_left(0.0, 0.0, 0.0);
    Point<3> upper_right(1.0, 1.0, 2.0);

    deallog << "Creating volume_mesh" << std::endl;
    GridGenerator::subdivided_hyper_rectangle(
        volume_mesh, n_subdivisions, bottom_left, upper_right);

    deallog << "Extracting boundary mesh" << std::endl;
    extract_surface(volume_mesh, surface_mesh);

    GridOut grid_out;
    GridOutFlags::Ucd ucd_flags;
    ucd_flags.write_preamble = true;
    grid_out.set_flags(ucd_flags);

    deallog << "Writing out surface_mesh to " << outfile << std::endl;
    std::ofstream out;
    out.open(outfile.c_str());
    grid_out.write_ucd(surface_mesh, out);
    out.close();

    return 0;
}

