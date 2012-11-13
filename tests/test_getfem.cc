#include <iostream>
#include <vector>
#include <typeinfo>

#include <getfem/bgeot_geometric_trans.h>
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_import.h>
#include <getfem/getfem_integration.h>
#include <getfem/getfem_mesh_fem.h>
#include <getfem/getfem_mesh_im.h>

#include <muParser.h>

#include <sloc/utils.h>

using namespace std;


template <class CONT>
void copy_to_cout(const CONT &c, const char *delim)
{
    std::copy(c.begin(), c.end(),
        std::ostream_iterator<typename CONT::value_type>(std::cout, delim));
}

double dot(bgeot::base_node x, bgeot::base_node y)
{
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

bgeot::base_node cross(bgeot::base_node x, bgeot::base_node y)
{
    //
    //  | i  j  k  |
    //  | x0 x1 x2 | 
    //  | y0 y1 y2 |
    //
    return bgeot::base_node(
                (x[1] * y[2] - y[1] * x[2]),
               -(x[0] * y[2] - y[0] * x[2]),
                (x[0] * y[1] - y[0] * x[1])
    );
}

bgeot::base_node calc_normal(bgeot::base_node A, bgeot::base_node B, bgeot::base_node C)
{
    bgeot::base_node n = cross(C-A, B-A);
    double mag = dot(n,n);
    return n / sqrt(mag);
}

bgeot::base_node calc_normal(bgeot::base_matrix G)
{
    // columns of G are the points of our triangle
    bgeot::base_node A(G(0,0), G(1,0), G(2,0));
    bgeot::base_node B(G(0,1), G(1,1), G(2,1));
    bgeot::base_node C(G(0,2), G(1,2), G(2,2));
    return calc_normal(A,B,C);
}

// ----------------------------------------------------------------------------

// XXX: conditionally compile this class in, only if <sloc/config.h> says that muparser is installed
class MuFunctor
{
public:
    MuFunctor(std::string fn)
    {
        p.DefineVar("x", &_x);
        p.DefineVar("y", &_y);
        p.DefineVar("z", &_z);
        p.SetExpr(fn.c_str());
    }

    double operator()(double x, double y, double z)
    {
        _x = x;
        _y = y;
        _z = z;
        return p.Eval();
    }

private:
    mu::Parser p;
    double _x, _y, _z;
};

// ----------------------------------------------------------------------------

void print_mesh_info(getfem::mesh& m)
{
    bgeot::base_node Pmin, Pmax;
    m.bounding_box(Pmin, Pmax);
    cout << "Mesh Bounding Box" << endl;
    cout << "  Pmin = " << Pmin << endl;
    cout << "  Pmax = " << Pmax << endl;
}

void print_mesh_points(getfem::mesh& m)
{
    for (dal::bv_visitor i(m.points_index()); !i.finished(); ++i)
    {
        cout << "Point of index " << i << ": " << m.points()[i] << endl;
    }
}

void print_mesh_elements(getfem::mesh& m)
{
    for (dal::bv_visitor cv(m.convex_index()); !cv.finished(); ++cv)
    {
        const int npts = m.structure_of_convex(cv)->nb_points();
        cout << "Element of index " << cv << " has " << npts << " points" << endl;

        if (true)
        {
            cout << "  points:\n\t";
            std::copy(m.points_of_convex(cv).begin(),
                      m.points_of_convex(cv).end(),
                      std::ostream_iterator<bgeot::base_node>(cout, "\n\t"));
            cout << "\r";
        }

        if (true)
        {
            cout << "  indices: ";
            std::copy(m.ind_points_of_convex(cv).begin(),
                      m.ind_points_of_convex(cv).end(),
                      std::ostream_iterator<bgeot::size_type>(cout, " "));
            cout << "\r\n";
        }

        if (true)
        {
            bgeot::pgeometric_trans pgt = m.trans_of_convex(cv);
            cout << "  geometric transformation: " << bgeot::name_of_geometric_trans(pgt) << endl;

            // convex of reference
            //bgeot::pconvex_ref cvr = pgt->convex_ref();

            //cout << "  reference normals:\n\t";
            //std::copy(cvr->normals().begin(), cvr->normals().end(), std::ostream_iterator<bgeot::base_small_vector>(cout, "\n\t"));
            //cout << "\r";
            //pgt->transform();

        }

        if (true)
        {
            bgeot::base_node pt = gmm::mean_value(m.points_of_convex(cv));
            cout << "  convex centroid: " << pt << endl;

            cout << "  convex area: " << m.convex_area_estimate(cv) << endl;

            //cout << "  convex quality: " << m.convex_quality_estimate(cv) << endl;
            //cout << "  convex radius: " << m.convex_radius_estimate(cv) << endl;
        }
    }
}

void print_integration_method()
{
    bgeot::pgeometric_trans pgt = bgeot::simplex_geotrans(2,1);
    cout << "Geometric Transformation: " << bgeot::name_of_geometric_trans(pgt) << endl;

    getfem::pintegration_method pim_exact = getfem::classical_exact_im(pgt);
    cout << "Exact Integration Method: " << getfem::name_of_int_method(pim_exact) << endl;

    getfem::pintegration_method pim_approx = getfem::classical_approx_im(pgt, 6);
    getfem::papprox_integration pai = pim_approx->approx_method();
    cout << "Approximate Integration Method: " << getfem::name_of_int_method(pim_approx) << endl;
    for (bgeot::size_type i = 0; i < pai->nb_points_on_convex(); ++i)
    {
        bgeot::scalar_type coeff = pai->coeff(i);
        bgeot::base_node point = pai->point(i);
        cout << "  " << coeff << " " << point << endl;
    }

}

void print_integral(getfem::mesh& mesh, std::string fn)
{

    // function we want to integrate!
    MuFunctor F(fn);

    bgeot::size_type N = mesh.dim();

    // finite element on mesh
    getfem::mesh_fem mf(mesh);
    mf.set_finite_element(getfem::fem_descriptor("FEM_PK(2,1)"));

    // integration method on mesh
    getfem::mesh_im mim(mesh);
    getfem::pintegration_method pim = getfem::classical_approx_im(bgeot::simplex_geotrans(2,1), 6);
    mim.set_integration_method(pim);
    
    // some useful vectors
    bgeot::base_matrix G;
    //std::vector<bgeot::scalar_type> coeff;


    // the final result
    double val = 0.0;


    // iterate through all the elements in mesh
    for (dal::bv_visitor cv(mesh.convex_index()); !cv.finished(); ++cv)
    {
        bgeot::pgeometric_trans pgt = mesh.trans_of_convex(cv);
        getfem::papprox_integration pai = getfem::get_approx_im_or_fail(mim.int_method_of_element(cv));
        getfem::pfem pf = mf.fem_of_element(cv);

        if (false)
        {
            // now, print out whatever is in mesh.points_of_convex(cv)
            //cout << typeid(mesh.points_of_convex(cv)).name() << endl;
            //cout << typeid(mesh.points_of_convex(cv).begin()).name() << endl;
            cout << "  mesh.points_of_convex(" << cv << "):\n\t";
            copy_to_cout(mesh.points_of_convex(cv), "\n\t");
            cout << "\r";

        }

        if (false)
        {
            cout << "  mf.nb_basic_dof_of_element(" << cv << ") = "
                 << mf.nb_basic_dof_of_element(cv) << endl;

            cout << "  mf.ind_basic_dof_of_element(" << cv << ") = ";
            copy_to_cout(mf.ind_basic_dof_of_element(cv), " ");
            cout << endl;
        }

        bgeot::vectors_to_base_matrix(G, mesh.points_of_convex(cv));
        //cout << "G = " << G << endl;

        bgeot::base_node normal = calc_normal(G);
        cout << "normal = " << normal << endl;

        //coeff.resize(mf.nb_basic_dof_of_element(cv));

        //gmm::copy(gmm::sub_vector(U, gmm::sub_index(mf.ind_basic_dof_of_element(cv))), coeff);

        getfem::fem_interpolation_context ctx(pgt, pf, bgeot::base_node(N), G, cv);

        for (unsigned int q = 0; q < pai->nb_points_on_convex(); ++q)
        {
            // set current point on reference convex (called xref)
            ctx.set_xref(pai->point(q));

            // evaluate shape functions at xref
            getfem::base_tensor t;
            ctx.base_value(t);

            // calculate x (geometric transform applied to xref)
            bgeot::base_node x;
            //x = pgt->transform(pai->point(q), mesh.points_of_convex(cv));
            x = pgt->transform(pai->point(q), G);

            if (true)
            {
                //cout << t.sizes() << endl;
                //cout << "value at " << pai->point(q) << " = " << t << endl;
                //cout << pai->point(q) << ";";
                //cout << "[" << t(0,0) << ", " << t(1,0) << ", " << t(2,0) << "];";
                //cout << x << endl;
                //cout << endl;
            }

            // update val
            double f = F(x[0], x[1], x[2]);
            val += f * pai->coeff(q) * ctx.J();
        }
    }

    cout << "Value of integral of " << fn << endl;
    cout << "  " << val << endl;

    cout << endl;
}

// ----------------------------------------------------------------------------

void test_load_native_mesh()
{
    string datadir = "/Users/luis/w/sloc/data";
    string meshfile = "small.mesh";
    string fullpath = datadir + string("/") + meshfile;

    getfem::mesh m;
    m.read_from_file(fullpath);

    //print_mesh_info(m);
    print_mesh_points(m);
    print_mesh_elements(m);

    //print_integral(m, "1");
    print_integral(m, "x*y");
}

void test_load_gid_mesh()
{
    string datadir = "/Users/luis/opt/getfem/getfem_toolbox/meshes";
    string meshfile = "tripod.GiD.msh";
    //string meshfile = "holed_disc_with_quadratic_2D_triangles.msh";
    string fullpath = datadir + string("/") + meshfile;
    string gidpath = string("gid:") + fullpath;
    string outfile = "test.mesh";

    getfem::mesh m;
    getfem::import_mesh(gidpath, m);

    print_mesh_info(m);
    print_mesh_points(m);
    print_mesh_elements(m);

    m.write_to_file(outfile);
    cout << "wrote " << outfile << endl;
}

// ----------------------------------------------------------------------------

int main(int argc, const char *argv[])
{
    cout << "----------------------------" << endl;
    //print_integration_method();
    test_load_native_mesh();
    //test_load_gid_mesh();
    return 0;
}
