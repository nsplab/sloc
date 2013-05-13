#include "sloc/bem_forward_problem_par.h"

using namespace std;
using namespace sloc;

// ----------------------------------------------------------------------------

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::valarray<T>& xs)
{
    const int n = xs.size();
    for (int i = 0; i < n; i++)
    {
        os << xs[i];
        if (i != n-1)
            os << " ";
    }
    return os;
}

// ----------------------------------------------------------------------------

BEM_Forward_Problem_P::BEM_Forward_Problem_P(const Parameters& parameters)
    : BEM_Forward_Problem(parameters)
{
    numprocs = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
}

BEM_Forward_Problem_P::~BEM_Forward_Problem_P()
{
}

void BEM_Forward_Problem_P::run()
{
    sloc::BEM_Forward_Problem::configure();

    // this class overrides the assemble_system() method, and builds
    // the system_matrix in parallel using multiple processes
    this->assemble_system();

    // process 0 is responsible for solving the linear system, and
    // and for saving the results to disk
    if (rank == 0)
    {
        sloc::BEM_Forward_Problem::solve_system();
        sloc::BEM_Forward_Problem::output_results();
    }
}

void BEM_Forward_Problem_P::assemble_range_contrib(unsigned int cv_begin, unsigned int cv_end, std::valarray<double>& contrib)
{
    using namespace dealii;

    const bgeot::size_type num_elts = surface_mesh.nb_convex();
    const unsigned int num_elts_by_proc = num_elts / numprocs;

    const unsigned int n_dofs = mf.nb_basic_dof();
    const unsigned int fe_dofs_per_cell = mf.nb_basic_dof_of_element(0);
    assert(fe_dofs_per_cell == 3);

    std::vector<unsigned int> local_dof_indices(fe_dofs_per_cell);
    dealii::Vector<double> local_matrix_row_i(fe_dofs_per_cell);

    // for holding the local points
    bgeot::base_matrix G;

    // loop indices
    unsigned int e, i, j, q;
    bgeot::size_type cv;

    // parallelized contribution from surface integral terms
    for (cv = cv_begin; cv < cv_end; ++cv)
    {
        e = cv - cv_begin;

        bgeot::pgeometric_trans pgt = surface_mesh.trans_of_convex(cv);
        getfem::papprox_integration pai = getfem::get_approx_im_or_fail(mim.int_method_of_element(cv));
        getfem::pfem pf = mf.fem_of_element(cv);

        bgeot::vectors_to_base_matrix(G, surface_mesh.points_of_convex(cv));

        // refactored here, since we know our elements are flat (revise this assumption later)
        dealii::Point<3> normal_q = triangle_normal(G);

        getfem::fem_interpolation_context ctx(pgt, pf, bgeot::base_node(3), G, cv);

        // get current element's dof indices
        //getfem::mesh_fem::ind_dof_ct elt_dof_indices = mf.ind_basic_dof_of_element(cv);

        // retrieve material data
        const unsigned int mat_id = material_data.get_material_id(cv);
        const double sigma_int = material_data.get_sigma_int(mat_id);
        const double sigma_ext = material_data.get_sigma_ext(mat_id);
        const double sigma_avg = (sigma_int + sigma_ext) / 2;
        const double K = (1.0 / (4 * numbers::PI)) * (sigma_int - sigma_ext) / sigma_avg;

        for (i = 0; i < n_dofs; ++i)
        {
            local_matrix_row_i = 0;

            bgeot::base_node p = mf.point_of_basic_dof(i);
            const Point<3> node_position(p[0], p[1], p[2]);

            for (q = 0; q < pai->nb_points_on_convex(); ++q)
            {
                // set current point on the reference elemnt
                ctx.set_xref(pai->point(q));

                // evaluate shape functions at xref (NOTE: this could be factored out?)
                bgeot::base_tensor t;
                ctx.base_value(t);

                // calculate x (geometric transform applied to xref)
                bgeot::base_node x;
                x = pgt->transform(pai->point(q), G);

                const Point<3> q_point(x[0], x[1], x[2]);
                const Point<3> R = q_point - node_position;
                const double R3 = std::pow(R.square(), 1.5);
                const double JxW = pai->coeff(q) * ctx.J();

                for (j = 0; j < fe_dofs_per_cell; ++j)
                {
                    const double shape_value_j = t(j,0);
                    const double term = K * (R / R3) * normal_q * shape_value_j * JxW;
                    local_matrix_row_i(j) += term;
                }
            }

            for (j = 0; j < fe_dofs_per_cell; ++j)
            {
                //system_matrix(i, elt_dof_indices[j]) += -local_matrix_row_i(j);
                contrib[j + fe_dofs_per_cell * (i + n_dofs * e)] = -local_matrix_row_i(j);
            }
        }
    }

    return;
}

void BEM_Forward_Problem_P::assemble_system_from_contrib(unsigned int cv_begin, unsigned int cv_end, std::valarray<double>& contrib, const int k)
{
    const unsigned int n_dofs = mf.nb_basic_dof();
    const unsigned int fe_dofs_per_cell = mf.nb_basic_dof_of_element(0);
    assert(fe_dofs_per_cell == 3);

    // loop indices
    unsigned int e, i, j, q;
    bgeot::size_type cv;

    for (cv = cv_begin; cv < cv_end; ++cv)
    {
        e = cv - cv_begin;
        getfem::mesh_fem::ind_dof_ct elt_dof_indices = mf.ind_basic_dof_of_element(cv);

        for (i = 0; i < n_dofs; ++i)
        {
            for (j = 0; j < fe_dofs_per_cell; ++j)
            {
                const double term = contrib[j + fe_dofs_per_cell * (i + n_dofs * e)];
                //cout << std::flush << "  k=" << k << " cv=" << cv << " i=" << i << " j=" << j << " M(" << i << "," << elt_dof_indices[j] << ") += " << term << endl;
                system_matrix(i, elt_dof_indices[j]) += term;
            }
        }
    }

    return;
}

void BEM_Forward_Problem_P::assemble_system()
{
    using namespace dealii;

    if (rank == 0)
        deallog << "BEM_Forward_Problem_P::assemble_system() T=" << timer.wall_time() << endl;

    bgeot::size_type num_elts = surface_mesh.nb_convex();
    unsigned int num_elts_by_proc = num_elts / numprocs;

    const unsigned int n_dofs = mf.nb_basic_dof();
    const int fe_dofs_per_cell = mf.nb_basic_dof_of_element(0);
    assert(fe_dofs_per_cell == 3);

    // loop indices
    unsigned int e, i, j, q;
    bgeot::size_type cv;

    // zero out system matrix, so we can call this method repeatedly
    system_matrix.reset_values();

    // print debug info, take care of system_rhs, and add identity to system_matrix
    if (rank == 0)
    {
        if (parameters.debug)
        {
            cout << std::flush;
            cout << "  num_elts / numprocs = " << num_elts_by_proc << endl;
            cout << "  num_elts % numprocs = " << (num_elts % numprocs) << endl;
            for (int k = 0; k < numprocs; ++k)
            {
                cout << "  range(" << k << ") = ["
                     << (k*num_elts_by_proc) << ","
                     << (k+1)*num_elts_by_proc << ")" << endl;
            }
            if (num_elts % numprocs != 0)
            {
                cout << "  range(" << numprocs << ") = ["
                     << numprocs * num_elts_by_proc << ","
                     << num_elts << "]" << endl;
            }
        }

        // Add up contribution from dipoles
        for (i = 0; i < n_dofs; ++i)
        {
            bgeot::base_node p = mf.point_of_basic_dof(i);
            dealii::Point<3> pt(p[0], p[1], p[2]);
            system_rhs(i) = 2 * dipole_sources.primary_contribution(pt);
        }

        // Add the identity matrix to system_matrix
        for (i = 0; i < n_dofs; ++i)
            system_matrix(i,i) += 1.0;
    }

    // array containing the integrals we'll need for assembling system_matrix
    std::valarray<double> kth_contrib;
    kth_contrib.resize(num_elts_by_proc * n_dofs * fe_dofs_per_cell);

    // break up the list of elements into a range
    unsigned int cv_begin = rank * num_elts_by_proc;
    unsigned int cv_end = (rank + 1) * num_elts_by_proc;

    // compute contributions on designated range of elements
    assemble_range_contrib(cv_begin, cv_end, kth_contrib);

    if (rank > 0)
    {
        // non-zero processes report results back: send message to process 0 (with tag 0)
        MPI::COMM_WORLD.Send(&kth_contrib[0], kth_contrib.size(), MPI::DOUBLE, 0, 0);
    }
    else if (rank == 0)
    {
        // what data do we currently own?
        //cout << std::flush << "owned by k=0\n\t" << kth_contrib << endl;

        // variable to denote the k-th process
        int k;

        //
        // Let's receive data from all other processes.
        //
        // Note that, in the first pass through this loop, we use the local
        // values of kth_contrib first, as computed by process 0. In the next
        // passes through the loop, we sequentially receive messages from the
        // other processes.
        //
        for (k = 0; k < numprocs; ++k)
        {
            if (k > 0)
            {
                // receive message from process k (with tag 0)
                MPI::COMM_WORLD.Recv(&kth_contrib[0], kth_contrib.size(), MPI::DOUBLE, k, 0);

                // what data did we receive?
                //cout << std::flush << "recv from k=" << k << "\n\t" << kth_contrib << endl;
            }

            // add k-th contribution to the system_matrix
            cv_begin = k * num_elts_by_proc;
            cv_end = (k + 1) * num_elts_by_proc;
            assemble_system_from_contrib(cv_begin, cv_end, kth_contrib, k);
        }

        // Handle remainder elements
        if (num_elts % numprocs != 0)
        {
            k = numprocs;
            cv_begin = numprocs * num_elts_by_proc;
            cv_end = num_elts;

            kth_contrib.resize((cv_end - cv_begin) * n_dofs * fe_dofs_per_cell);
            assemble_range_contrib(cv_begin, cv_end, kth_contrib);
            assemble_system_from_contrib(cv_begin, cv_end, kth_contrib, k);
        }
    }

    if ((rank == 0) && debug)
    {
        sloc::write_vector("system_rhs.dat", system_rhs);
        sloc::write_matrix("system_matrix.dat", system_matrix);
    }

    return;
}

