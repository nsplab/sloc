
#include <sloc/utils.h>
#include <getfem/bgeot_kdtree.h>
#include <getfem/dal_bit_vector.h>

using std::cout;
using std::cerr;
using std::endl;

using bgeot::base_node;
using bgeot::size_type;
using bgeot::dim_type;

bool quick = false;

void brute_force_points_in_box(const std::vector<base_node>& v,
                               dal::bit_vector& ipts,
                               const base_node& bmin,
                               const base_node& bmax)
{
    ipts.clear();
    for (size_type i = 0; i < v.size(); ++i) {
        bool is_in = true;
        for (size_type k = 0; k < bmin.size(); ++k) {
            if (v[i].at(k) < bmin[k] || v[i].at(k) > bmax[k]) {
                is_in = false;
                break;
            }
        }
        if (is_in)
            ipts.add(i);
    }
}

void verify_points_in_box(const std::vector<base_node>& v,
                          bgeot::kdtree& tree,
                          const base_node& bmin, 
                          const base_node& bmax)
{
    //cout << "verify_points_in_box: " << bmin << ", " << bmax << endl;

    dal::bit_vector bv1, bv2;

    brute_force_points_in_box(v, bv1, bmin, bmax);

    bgeot::kdtree_tab_type ipts;
    tree.points_in_box(ipts, bmin, bmax);

    for (size_type i = 0; i < ipts.size(); ++i)
        bv2.add(ipts[i].i);

    if (bv1 != bv2) {
        cout << "verify_points_in_box: error, brute force gave points\n" << bv1
             << ",\nwhile points_in_box returned: " << bv2 << endl;
    } else {
        cout << ".";
        cout.flush();
    }
    assert(bv1 == bv2);
}

void check_tree()
{
    bgeot::kdtree tree;
    std::vector<size_type> ipts;
    std::vector<base_node> pts;

    base_node one(2); gmm::fill(one, 1);

    pts.push_back(base_node(0.4, 0.7));
    pts.push_back(base_node(1.4, 0.8));
    pts.push_back(base_node(2.4, 0.3));
    pts.push_back(base_node(3.4, 0.1));
    pts.push_back(base_node(0.5, 2.2));
    pts.push_back(base_node(-0.8, 0.4));

    for (size_type i = 0; i < 50; ++i)
        pts.push_back(base_node(gmm::random(), gmm::random()));

    for (size_type i = 0; i < pts.size(); ++i)
        tree.add_point(pts[i]);

    /* all points */
    verify_points_in_box(pts, tree, base_node(0.0, 0.0), base_node(6.0, 6.0));

    /* point search */
    for (size_type i = 0; i < pts.size(); ++i) {
        for (size_type j = 0; j < 2; ++j) {
            base_node bmin = pts[i], bmax = pts[i];
            if (j == 0) {
                bmin -= 1e-5 * one;
                bmax += 1e-5 * one;
            }
            verify_points_in_box(pts, tree, bmin, bmax);
        }
    }

    /* check special cases (all points the same, all points on same plane etc) */
    for (size_type repeat = 0; repeat < 200; ++repeat) {
        pts.clear();
        tree.clear();

        double pval[3] = { 0.3, 24.0, -2.3 };
        bool sameplane[3];

        for (size_type i = 0; i < 3; ++i)
            sameplane[i] = ((rand() % 3) == 0);

        for (size_type i = 0; i < 200; ++i) {
            base_node pt(3);
            for (size_type k = 0; k < pt.size(); ++k) {
                if (sameplane[k])
                    pt[k] = pval[1];
                else
                    pt[k] = pval[rand() % 3];
            }
            pts.push_back(pt);
            tree.add_point(pt);
        }

        verify_points_in_box(pts, tree, base_node(-100.0, -100.0, -100.0), base_node(100.0, 100.0, 100.0));
        verify_points_in_box(pts, tree, base_node(-100.0, -100.0, -100.0), base_node(0.0, 0.0, 0.0));
        verify_points_in_box(pts, tree, pts[0], pts[0]);
        verify_points_in_box(pts, tree, pts[0], pts[1]);
    }
    cout << endl << "the kdtree is ok!" << endl;
}

void speed_test(unsigned int N, unsigned int NPT, unsigned int nrepeat)
{
    bgeot::kdtree tree;
    base_node pt(N);

    cout << "speed test for the kdtree\n";
    double t = gmm::uclock_sec();
    for (size_type i = 0; i < NPT; ++i) {
        for (dim_type k = 0; k < N; ++k)
            pt[k] = gmm::random(double()) * 2.0;
        tree.add_point(pt);
        assert(pt.refcnt() > 1);
    }
    cout << "point list built in " << gmm::uclock_sec() - t << " seconds.\n";

    t = gmm::uclock_sec();
    bgeot::kdtree_tab_type ipts;
    tree.points_in_box(ipts, pt, pt);
    cout << "tree built in " << gmm::uclock_sec() - t << " seconds.\n";
    
    t = gmm::uclock_sec();
    base_node bmin(0.25, 0.25, 0.4), bmax(0.5, 0.5, 3.3);
    size_type npt = 0;
    for (size_type c = 0; c < nrepeat; ++c) {
        tree.points_in_box(ipts, bmin, bmax);
        npt = ipts.size();
    }
    cout << "BOX QUERY: nb points in " << bmin << ":" << bmax << " is : " << npt << endl;
    cout << "average query time is: " << (gmm::uclock_sec()-t)/nrepeat*1e6 << " microseconds\n";
}

int main(int argc, const char *argv[])
{
    if (argc == 2 && strcmp(argv[1], "--quick") == 0)
        quick = true;

    check_tree();

    if (!quick)
        speed_test(3, 300000, 20000);
    else
        speed_test(2, 10000, 100);

    return 0;
}
