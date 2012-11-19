#include "sloc/bem_forward_problem_par.h"

using namespace std;
using namespace sloc;

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

    // this class overrides the assemble_system() method, to build
    // the system_matrix in parallel
    this->assemble_system();

    // process 0 is responsible for solving the linear system, and
    // and for saving the results to disk
    if (rank == 0)
    {
        if (true)
        {
            sloc::write_vector("system_rhs.dat", system_rhs);
            sloc::write_matrix("system_matrix.dat", system_matrix);
        }

        //sloc::BEM_Forward_Problem::solve_system();
        //sloc::BEM_Forward_Problem::output_results();
    }
}

void BEM_Forward_Problem_P::assemble_system()
{
    using namespace std;
    using namespace dealii;

    if (rank == 0)
        deallog << "BEM_Forward_Problem_P::assemble_system() T=" << timer.wall_time() << endl;

    const unsigned int n_dofs = mf.nb_basic_dof();
    const int fe_dofs_per_cell = mf.nb_basic_dof_of_element(0);
    assert(fe_dofs_per_cell == 3);

    std::vector<unsigned int> local_dof_indices(fe_dofs_per_cell);
    dealii::Vector<double> local_matrix_row_i(fe_dofs_per_cell);

    std::valarray<double> kth_contrib;
    kth_contrib.resize(n_dofs * fe_dofs_per_cell);
    kth_contrib = 0;

    unsigned int i, j, q;

    // for holding the local points
    bgeot::base_matrix G;

    // zero out system matrix, so we can call this method repeatedly
    system_matrix.reset_values();

    bgeot::size_type num_elts = surface_mesh.nb_convex();
    unsigned int num_elts_by_proc = num_elts / numprocs;
    unsigned int cv_begin = rank * num_elts_by_proc;
    unsigned int cv_end = (rank + 1) * num_elts_by_proc;
    bgeot::size_type cv;

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

    // parallelized contribution from surface integral terms
    for (cv = cv_begin; cv < cv_end; ++cv)
    {
        bgeot::pgeometric_trans pgt = surface_mesh.trans_of_convex(cv);
        getfem::papprox_integration pai = getfem::get_approx_im_or_fail(mim.int_method_of_element(cv));
        getfem::pfem pf = mf.fem_of_element(cv);

        bgeot::vectors_to_base_matrix(G, surface_mesh.points_of_convex(cv));

        // refactored here, since we know our elements are flat (revise this assumption later)
        dealii::Point<3> normal_q = triangle_normal(G);

        getfem::fem_interpolation_context ctx(pgt, pf, bgeot::base_node(3), G, cv);

        // get current element's dof indices
        getfem::mesh_fem::ind_dof_ct elt_dof_indices = mf.ind_basic_dof_of_element(cv);

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
                kth_contrib[i * fe_dofs_per_cell + j] = -local_matrix_row_i(j);
                //kth_contrib[i * fe_dofs_per_cell + j] = rank;
            }
        }
    }

    // Report back to process 0
    if (rank > 0)
    {
        MPI::COMM_WORLD.Send(&kth_contrib[0], kth_contrib.size(), MPI::DOUBLE, 0, 0);
    }
    else if (rank == 0)
    {

        // own data
        cout << "own k=0: ";
        for (i = 0; i < kth_contrib.size(); ++i) {cout << kth_contrib[i] << " ";} cout << endl;

        // Receive data from all other processes
        // (note that own contribution is added directly to system_matrix first)
        for (int k = 0; k < numprocs; ++k)
        {
            // receive k-th contribution
            if (k > 0)
            {
                kth_contrib = 0;

                // mpi receive
                MPI::COMM_WORLD.Recv(&kth_contrib[0], kth_contrib.size(), MPI::DOUBLE, k, 0);
                if (true)
                {
                    //cout << "Received data from process " << k << "\n\t";
                    cout << "recv k=" << k << ": ";
                    for (i = 0; i < kth_contrib.size(); ++i) {cout << kth_contrib[i] << " ";} cout << endl;
                }
            }

            // add k-th contribution to the system_matrix
            cv_begin = k * num_elts_by_proc;
            cv_end = (k + 1) * num_elts_by_proc;
            for (cv = cv_begin; cv < cv_end; ++cv)
            {
                getfem::mesh_fem::ind_dof_ct elt_dof_indices = mf.ind_basic_dof_of_element(cv);
                for (i = 0; i < n_dofs; ++i) {
                    for (j = 0; j < fe_dofs_per_cell; ++j) {
                        const double term = kth_contrib[i * fe_dofs_per_cell + j];
                        cout << "  k=" << k << " cv=" << cv << " i=" << i << " j=" << j << " M(" << i << "," << elt_dof_indices[j] << ") += " << term << endl;
                        system_matrix(i, elt_dof_indices[j]) += term;
                        //system_matrix(i, elt_dof_indices[j]) += k * std::pow(10., j*1.);
                    }
                }
            }
        }

        // Handle remainder elements
        cv_begin = numprocs * num_elts_by_proc;
        cv_end = num_elts;
        for (cv = cv_begin; cv < cv_end; ++cv)
        {
            getfem::mesh_fem::ind_dof_ct elt_dof_indices = mf.ind_basic_dof_of_element(cv);
            for (i = 0; i < n_dofs; ++i) {
                for (j = 0; j < fe_dofs_per_cell; ++j) {
                    const double term = kth_contrib[i * fe_dofs_per_cell + j];
                    system_matrix(i, elt_dof_indices[j]) += term;
                }
            }
        }
    }

    return;
}

