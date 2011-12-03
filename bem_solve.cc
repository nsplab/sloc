
/* deal.II includes */

// for handling triangulations
#include <deal.II/grid/tria.h>

// for looping over cells and faces
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>

// input of grids in various formats
#include <deal.II/grid/grid_in.h>

// output of grids in various graphics formats
#include <deal.II/grid/grid_out.h>

// for GridGenerator::hyper_rectangle function
#include <deal.II/grid/grid_generator.h>

// various predefined boundaries (HyperBallBoundary, ...)
#include <deal.II/grid/tria_boundary_lib.h>

// linear algebra includes
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

// for managing the distribution and number of the degrees of freedom
// (i.e., associating dofs with vertices, lines, and cells).
#include <deal.II/dofs/dof_handler.h>

// provide information about the degrees of freedom local to a cell
#include <deal.II/dofs/dof_accessor.h>

// for the function DoFTools::map_dofs_to_support_points
#include <deal.II/dofs/dof_tools.h>

// for the Lagrange finite elements we'll use
#include <deal.II/fe/fe_q.h>

// for the Gauss-Legendre quadrature rule we'll be using
#include <deal.II/base/quadrature_lib.h>

// for evaluating the finite element shape functions at quadrature points of cell
#include <deal.II/fe/fe_values.h>

// for mapping reference element to arbitrary cell
#include <deal.II/fe/mapping_q.h>

// for writing out the output file given the shape function coefficients from our solution
#include <deal.II/numerics/data_out.h>

// for the function VectorTools::integrate_difference
#include <deal.II/numerics/vectors.h>

// various utility functions
#include <deal.II/base/utilities.h>
#include <deal.II/base/timer.h>

#include <deal.II/base/table.h>

/* standard includes */
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

/* other functions */
#include "misc_io.h"

/* namespaces */
using namespace dealii;


// ----------------------------------------------------------------------------

struct DipoleSource
{
    double primary_contribution(const Point<3> &field_position);
    Point<3> position;
    Point<3> dipole;
};

double DipoleSource::primary_contribution(const Point<3> &field_position)
{
    const double sigma0 = 1.0;
    Point<3> R = field_position - position;
    return (1.0/(4 * numbers::PI * sigma0)) * (dipole * R) / std::pow(R.square(), 1.5);
}

// ----------------------------------------------------------------------------

class MaterialData
{
    public:
        MaterialData(const unsigned int n_layers);

        double get_sigma_int(const unsigned int material_id) const;

        double get_sigma_ext(const unsigned int material_id) const;

        double get_sigma_avg(const unsigned int material_id) const;

  //void _init(const double data[][2]);

        void _init(const double data[][2], const double average_data[][1]);

    private:
        const unsigned int n_layers;
        Table<2,double> sigma;
        Table<2,double> average_sigma;
};

//Conductivities stored in a table, interior conductivity in first column.
MaterialData::MaterialData(const unsigned int n_layers)
    :
    n_layers(n_layers),
    sigma(n_layers, 2),
    average_sigma(n_layers, 2)
{}

double MaterialData::get_sigma_int(const unsigned int material_id) const
{
    Assert(material_id < n_layers,
           ExcIndexRange(material_id, 0, n_layers));
    return sigma[material_id][0];
}

double MaterialData::get_sigma_ext(const unsigned int material_id) const
{
    Assert(material_id < n_layers,
           ExcIndexRange(material_id, 0, n_layers));
    return sigma[material_id][1];
}

//void MaterialData::_init(const double data[][2])
//{
//    for(unsigned int m = 0; m < n_layers; ++m)
//    {
//        sigma[m][0] = data[m][0];
//        sigma[m][1] = data[m][1];
//    }
//}
double MaterialData::get_sigma_avg(const unsigned int material_id) const
{
    Assert(material_id < n_layers,
           ExcIndexRange(material_id, 0, n_layers));
    return average_sigma[material_id][0];
}

void MaterialData::_init(const double data[][2], const double average_data[][1])
{
    for(unsigned int m = 0; m < n_layers; ++m)
    {
        sigma[m][0] = data[m][0];
        sigma[m][1] = data[m][1];
        average_sigma[m][0] = average_data[m][0];
    }
}

//-----------------------------------------------------------------------

class BEM_ForwardProblem
{
public:
    BEM_ForwardProblem();
    void run();
    ~BEM_ForwardProblem();

private:

    void configure();
    void compute_area();
    void assemble_system();
    void solve_system();
    void output_results();
    void compute_exterior_solution();

    MaterialData *material_data;

    // Location and components of dipole source
    DipoleSource dipole_source;
    static const double sigma0 = 1.0;
  //static const double sigma_int = 200.0;
  //static const double sigma_ext = 0.0;
  //double sigma_avg;

    // triangulation of domain: use 2d cells embedded in 3d space
    Triangulation<2,3> tria;

    // use a 2d Lagrange finite element embedded in 3d space
    FE_Q<2,3> fe;

    // the degrees of freedom for 2d elements embedded in 3d space
    DoFHandler<2,3> dh;

    // the map function between 2d reference element to its cell in 3d space 
    MappingQ<2,3> mapping;
    
    // the 2d quadrature rule to use when integrating our 2d reference element
    Quadrature<2> quadrature;

    // The system matrix (I-C), and the rhs vector g.
    FullMatrix<double> system_matrix;
    Vector<double> system_rhs;

    // The solution $\phi$
    Vector<double> phi;

    // parameters for the controlling linear solver (tolerance, logging, etc.)
    SolverControl solver_control;

    // wallclock time
    Timer timer;
};

BEM_ForwardProblem::BEM_ForwardProblem()
    :
    material_data(0),
    fe(1),                      // fe degree 1
    dh(tria),                   // attach triangulation to our dof_handler
    mapping(1, true)           // mapping degree 1 (also, use same mapping on interior cells)
{
    timer.start();
}

BEM_ForwardProblem::~BEM_ForwardProblem()
{
    if(material_data)
    {
      delete material_data;
    }
}

void BEM_ForwardProblem::configure()
{
    using namespace std;

    deallog << "BEM_ForwardProblem::configure() " << timer.wall_time() << endl;

    // configure the dipole source
    const double d = pow(1/3.0, 0.5);
    dipole_source.position = Point<3>(0.5, 0.5, 0.5);
    dipole_source.dipole = Point<3>(d, d, d);
    //sigma_avg = (sigma_int + sigma_ext) / 2.0;

    /*
    // Input array from text file. For now, adjust the number of lines
    // and the file name here.
    lines = 3
    ifstream input_file("sigmafile.txt");
    double temp_sigma[lines][2]
    if(input_file.is_open())
    {
      for(int i = 0; !input_file.eof(); ++i)
        {
          input_file >> temp_sigma[i][0];
          input_file >> temp_sigma[i][1];
        }
      input_file.close();
    }
    else {
        cerr << "Conductivity file not found." < endl;
        exit(1);
    }
    const double sigma_data = temp_sigma
    */

    // Sample array of conductivities with 3 surfaces.
    const double sigma_data[3][2] = {
      {0, 2000},
      {2000, 10},
      {10, 0}
    }; 

    const double avg_sigma_data[3][1] = {
      {0},
      {1.5},
      {100},
    };

    material_data = new MaterialData(3);
    //    material_data->_init(sigma_data);
    material_data->_init(sigma_data, avg_sigma_data);

    // use Gauss-Legendre quadrature rule of order 4
    quadrature = QGauss<2>(4);

    // configure the solver control
    solver_control.log_frequency(1);
    solver_control.log_history(false);
    solver_control.log_result(true);
    solver_control.set_max_steps(100);
    solver_control.set_tolerance(1.0e-10);

    // read the boundary mesh for our domain
    read_ucd_mesh("doublespheresurf.ucd", tria);
    write_triangulation("sphere.inp", tria);

    // initialize vectors
    dh.distribute_dofs(fe);
    const unsigned int n_dofs = dh.n_dofs();
    system_matrix.reinit(n_dofs, n_dofs);
    system_rhs.reinit(n_dofs);
    phi.reinit(n_dofs);
}

void BEM_ForwardProblem::assemble_system()
{
    deallog << "BEM_ForwardProblem::assemble_system() " << timer.wall_time() << std::endl;

    FEValues<2,3> fe_v(mapping, fe, quadrature,
        update_values |
        update_cell_normal_vectors |
        update_quadrature_points |
        update_JxW_values);

    const unsigned int n_dofs = dh.n_dofs();
    const unsigned int n_q_points = fe_v.n_quadrature_points;

    std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);

    Vector<double> local_matrix_row_i(fe.dofs_per_cell);

    std::vector<Point<3> > support_points(n_dofs);
    DoFTools::map_dofs_to_support_points<2,3>(mapping, dh, support_points);
    //write_points("support_points.dat", support_points);

    DoFHandler<2,3>::active_cell_iterator
        cell = dh.begin_active(),
        endc = dh.end();

    double sigma_avg;

    for (unsigned int i = 0; i < n_dofs; ++i)
    {
        system_rhs(i) = 2 * dipole_source.primary_contribution(support_points[i]);
    }

    //std::ofstream q_out, n_out;
    //q_out.open("qpts.dat");
    //n_out.open("normals.dat");

    for (cell = dh.begin_active(); cell != endc; ++cell)
    {
        fe_v.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        const std::vector<Point<3> > &q_points = fe_v.get_quadrature_points();
        const std::vector<Point<3> > &normals = fe_v.get_cell_normal_vectors();

        //write_points(q_out, q_points);
        //write_points(n_out, normals);

        const double sigma_int = material_data->get_sigma_int(cell->material_id());
        const double sigma_ext = material_data->get_sigma_ext(cell->material_id());
        sigma_avg = (sigma_int + sigma_ext)/2;

        for (unsigned int i = 0; i < n_dofs; ++i)
        {
            local_matrix_row_i = 0;

            for (unsigned int q = 0; q < n_q_points; ++q)
            {
                const double K = (1.0/(4 * numbers::PI)) * ((sigma_int - sigma_ext) / sigma_avg);
                const Point<3> R = q_points[q] - dipole_source.position;
                double R3 = std::pow(R.square(), 1.5);

                for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
                    local_matrix_row_i(j) += K * (R / R3) * normals[q] * fe_v.shape_value(j,q) * fe_v.JxW(q);
            }

            for (unsigned int j = 0; j < fe.dofs_per_cell; ++j)
                system_matrix(i,local_dof_indices[j]) += -local_matrix_row_i(j);
        }
    }

    //q_out.close();
    //n_out.close();

    for (unsigned int i = 0; i < n_dofs; ++i)
        system_matrix(i,i) += 1.0;

    //write_matrix("system_matrix.dat", system_matrix);
    //write_vector("system_rhs.dat", system_rhs);
}

void BEM_ForwardProblem::solve_system()
{
    deallog << "BEM_ForwardProblem::solve_system() " << timer.wall_time() << std::endl;
    SolverGMRES<Vector<double> > solver(solver_control);
    solver.solve(system_matrix, phi, system_rhs, PreconditionIdentity());
}

void BEM_ForwardProblem::output_results()
{
    deallog << "BEM_ForwardProblem::output_results() " << timer.wall_time() << std::endl;

    DataOut<2, DoFHandler<2,3> > data_out;
    data_out.attach_dof_handler(dh);
    data_out.add_data_vector(phi, "phi");
    data_out.build_patches(mapping, mapping.get_degree());

    std::ofstream file("phi.vtk");
    data_out.write_vtk(file);
}


void BEM_ForwardProblem::compute_exterior_solution()
{
    deallog << "BEM_ForwardProblem::compute_exterior_solution() " << timer.wall_time() << std::endl;

    Triangulation<3> external_tria;
    read_ucd_mesh("doublesphere.ucd", external_tria);

    FE_Q<3>         external_fe(1);
    DoFHandler<3>   external_dh(external_tria);
    Vector<double>  external_phi;

    //    external_tria.refine_global(5);
    external_dh.distribute_dofs(external_fe);
    external_phi.reinit(external_dh.n_dofs());

    DoFHandler<2,3>::active_cell_iterator
        cell = dh.begin_active(),
        endc = dh.end();

    FEValues<2,3> fe_v(mapping, fe, quadrature,
        update_values |
        update_cell_normal_vectors |
        update_quadrature_points |
        update_JxW_values);

    const unsigned int n_q_points = fe_v.n_quadrature_points;

    std::vector<unsigned int> dofs(fe.dofs_per_cell);

    std::vector<double> local_phi(n_q_points);
    std::vector<Point<3> > external_support_points(external_dh.n_dofs());
    DoFTools::map_dofs_to_support_points<3>(StaticMappingQ1<3>::mapping, external_dh, external_support_points);

    double sigma_avg;

    for (unsigned int i = 0; i < external_dh.n_dofs(); ++i)
    {
        external_phi(i) = 2 * dipole_source.primary_contribution(external_support_points[i]);
    }

    for (cell = dh.begin_active(); cell != endc; ++cell)
    {
        fe_v.reinit(cell);

        const std::vector<Point<3> > &q_points = fe_v.get_quadrature_points();
        const std::vector<Point<3> > &normals = fe_v.get_cell_normal_vectors();

        const double sigma_int = material_data->get_sigma_int(cell->material_id());
        const double sigma_ext = material_data->get_sigma_ext(cell->material_id());
        sigma_avg = (sigma_int + sigma_ext) / 2;

        cell->get_dof_indices(dofs);
        fe_v.get_function_values(phi, local_phi);

        for (unsigned int q = 0; q < n_q_points; ++q)
        {
            const double K = (1.0/(4 * numbers::PI)) * ((sigma_int - sigma_ext) / sigma_avg);
            const Point<3> R = q_points[q] - dipole_source.position;
            double R3 = std::pow(R.square(), 1.5);

            for (unsigned int i = 0; i < external_dh.n_dofs(); ++i)
                external_phi(i) += K * (local_phi[q] * (R / R3) * normals[q]) * fe_v.JxW(q);
        }
    }
    // Timer after loops.
    //    deallog << "Writing external image file: " << timer.wall_time() << std::endl;
    // Number of degrees of freedom:
    //deallog << "Number of degrees of freedom:" << external_dh.n_dofs() << std::endl;
    // Number of quadrature points:
    //deallog << "Number of q_points:" << n_q_points << std::endl;


    std::ofstream file("external_phi.vtk");
    DataOut<3> data_out;
    data_out.attach_dof_handler(external_dh);
    data_out.add_data_vector(external_phi, "external_phi");
    data_out.build_patches();
    data_out.write_vtk(file);
}


// just a test.
// we calculate the area by integrating the scalar 1
// over the entire triangulation.
void BEM_ForwardProblem::compute_area()
{
    deallog << "BEM_ForwardProblem::compute_area() " << timer.wall_time() << std::endl;

    double area = 0;
    const double one = 1;

    FEValues<2,3> fe_v(mapping, fe, quadrature, update_JxW_values);

    const unsigned int n_q_points = fe_v.n_quadrature_points;

    std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);

    DoFHandler<2,3>::active_cell_iterator
        cell = dh.begin_active(),
        endc = dh.end();

    for (cell = dh.begin_active(); cell != endc; ++cell)
    {
        fe_v.reinit(cell);
        cell->get_dof_indices(local_dof_indices);

        for (unsigned int q = 0; q < n_q_points; ++q)
            area += one * fe_v.JxW(q);
    }

    deallog << "triangulation area = " << area << std::endl;
}

void BEM_ForwardProblem::run()
{
    deallog << "BEM_ForwardProblem::run() " << timer.wall_time() << std::endl;
    configure();
    compute_area();
    assemble_system();
    solve_system();
    output_results();
    compute_exterior_solution();
    deallog << "last timer = " << timer.wall_time() << std::endl;
}

int main(void)
{
    deallog.depth_console(3);

    deallog << "main()" << std::endl;

    try
    {
        BEM_ForwardProblem fwd_problem;
        fwd_problem.run();
    }
    catch (std::exception &exc)
    {
        using namespace std;
        string sep = "----------------------------------------------------";
        cerr << endl
             << sep << endl
             << "Exception on processing: " << endl
             << exc.what() << endl
             << "Aborting!" << endl
             << sep << endl;
        return 1;
    }
    catch (...)
    {
        using namespace std;
        string sep = "----------------------------------------------------";
        cerr << endl
             << sep << endl
             << "Unknown exception!" << endl
             << "Aborting!" << endl
             << sep << endl;
        return 1;
    }

    return 0;
}
