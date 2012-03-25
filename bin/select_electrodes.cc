/*
 * Given a mesh and a potentials data file, randomly select electrodes from the
 * desired layers and write out an electrodes file we can use in bem_inverse_solve.
 */

#include "mesh.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/smart_ptr.hpp>

using namespace std;
using namespace boost::algorithm;

typedef unsigned int id;
typedef map<id,id> map_layers_t;

typedef set<id> set_t;
typedef boost::shared_ptr<set_t> shared_set_t;
typedef map<id,shared_set_t> map_id_set_t;

typedef vector<id> vector_t;
typedef boost::shared_ptr<vector_t> shared_vector_t;
typedef map<id,shared_vector_t> map_id_vector_t;

// ----------------------------------------------------------------------------

bool split_layer_num(std::string str, unsigned int &layer, unsigned int &num)
{
    vector<string> pair;
    split(pair, str, is_any_of("/:"), token_compress_on);

    if (pair.size() != 2)
        return false;

    try {
        layer = boost::lexical_cast<unsigned int>(pair[0]);
        num = boost::lexical_cast<unsigned int>(pair[1]);
    } catch (boost::bad_lexical_cast &) {
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------

void usage(const char *pgm)
{
    cerr << "Usage: " << pgm
         << " -m meshfile -p datfile -o outfile layer1/count1 layer2/count2 ..."
         << endl;
    exit(1);
}

void process_args(int argc, char *argv[], string& meshfile, string& datfile, string& outfile, map_layers_t& layers)
{
    char *pgm = argv[0];

    if (argc == 1)
        usage(pgm);

    meshfile = "";
    datfile = "";
    outfile = "";
    layers.clear();

    while (argc > 1)
    {
        if (strncmp(argv[1], "-m", 3) == 0)
        {
            // read the next argument (expecting meshfile)
            argc--; argv++;
            if (argc > 1)
                meshfile = argv[1];
            else
                break;
        }
        else if (strncmp(argv[1], "-p", 3) == 0)
        {
            // read the next argument (expecting datfile)
            argc--; argv++;
            if (argc > 1)
                datfile = argv[1];
            else
                break;
        }
        else if (strncmp(argv[1], "-o", 3) == 0)
        {
            // read the next argument (expecting outfile)
            argc--; argv++;
            if (argc > 1)
                outfile = argv[1];
            else
                break;
        }
        else
        {
            // assume it's a layer string
            unsigned int mat, num;
            bool valid = split_layer_num(argv[1], mat, num);

            if (!valid)
            {
                cerr << "Invalid layer count '" << argv[1] << "'" << endl;
                usage(pgm);
            }

            if (layers.count(mat) > 0)
            {
                cerr << "Duplicate layer specified: '" << argv[1] << "'" << endl;
                usage(pgm);
            }

            layers[mat] = num;
        }

        // advance to the next argument
        argc--; argv++;
    }

    if (meshfile.empty())
    {
        cerr << "Missing mesh file!" << endl;
        usage(pgm);
    }

    if (datfile.empty())
    {
        cerr << "Mssing data file!" << endl;
        usage(pgm);
    }

    if (outfile.empty())
    {
        cerr << "Missing output file!" << endl;
        usage(pgm);
    }

    if (layers.empty())
    {
        cerr << "No layers selected!" << endl;
        usage(pgm);
    }
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    int i,j;
    unsigned int mat, num;
    unsigned int total_num;

    // process cli args
    string meshfile, datfile, outfile;
    map_layers_t layers;
    process_args(argc, argv, meshfile, datfile, outfile, layers);

    // read the mesh
    sloc::Mesh mesh;
    cout << "Reading mesh '" << meshfile << "'" << endl;
    mesh.read_ucd(meshfile.c_str());

    // read the data file
    vector<double> phi;
    cout << "Reading potentials '" << datfile << "'" << endl;
    ifstream in;
    in.open(datfile.c_str());
    if (!in.is_open())
    {
        cerr << "Could not open potentials data file '" << datfile << "'" << endl;
        usage(argv[0]);
    }
    while (!in.eof())
    {
        string line;

        getline(in, line);
        if (line.empty())
            break;

        try {
            phi.push_back(boost::lexical_cast<double>(line));
        } catch (boost::bad_lexical_cast &) {
            cerr << "Bad line in potentials data file: " << line << endl;
            exit(2);
        }
    }
    in.close();

    // only cells are associated with materials, but we want an association
    // between nodes and material ids. so, we need to walk through each of the
    // cell's nodes and build a map from material ids to its set of nodes
    map_id_set_t nodes_by_mat;
    long cell[mesh.n_cell_nodes()];
    for (i = 0; i < mesh.n_cells(); i++)
    {
        int cell_mat_id;
        mesh.get_mat(i, cell_mat_id);

        if (layers.count(cell_mat_id) > 0)
        {
            mesh.get_cell(i, cell);

            if (nodes_by_mat.count(cell_mat_id) == 0)
                nodes_by_mat[cell_mat_id] = shared_set_t(new set_t);

            for (j = 0; j < mesh.n_cell_nodes(); j++)
                nodes_by_mat[cell_mat_id]->insert(cell[j]);
        }
    }

    // add up the layer counts
    // (1) only add counts whose layers that are present in the mesh file
    // (2) if desired count is larger than the number of actual nodes in the mesh
    //     by that material id, then use the actual number as the cutoff value
    total_num = 0;
    for (map_layers_t::iterator it = layers.begin(); it != layers.end(); ++it)
    {
        mat = it->first;
        num = it->second;
        if (nodes_by_mat.count(mat) > 0)
        {
            if (num > nodes_by_mat[mat]->size())
                it->second = nodes_by_mat[mat]->size();
            total_num += it->second;
        }
    }

    // now that duplicates are removed, let's copy the sets of nodes
    // into corresponding vectors, so we can shuffle them later
    map_id_vector_t nodes;
    for (map_id_set_t::iterator it = nodes_by_mat.begin(); it != nodes_by_mat.end(); ++it)
    {
        mat = it->first;

        if (nodes.count(mat) == 0)
            nodes[mat] = shared_vector_t(new vector_t);

        for (set_t::iterator sit = it->second->begin(); sit != it->second->end(); ++sit)
            nodes[mat]->push_back(*sit);
    }
    nodes_by_mat.clear();

    // shuffle each of the vectors of nodes and copy the first num elements
    srand(time(NULL));
    map_id_vector_t random_nodes;
    for (map_id_vector_t::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        unsigned int mat = it->first;
        unsigned int num = layers[mat];

        //cout << "shuffling layer " << mat << endl;
        random_shuffle(it->second->begin(), it->second->end());

        if (random_nodes.count(mat) == 0)
            random_nodes[mat] = shared_vector_t(new vector_t);

        // copy the first num entries in shuffled vector
        cout << "Taking " << num << " random nodes from layer " << mat << endl;
        unsigned count = 0;
        for (vector_t::iterator vit = it->second->begin(); vit != it->second->end(); ++vit)
        {
            random_nodes[mat]->push_back(*vit);
            if (++count >= num) break;
        }
    }
    nodes.clear();

    // write out the electrodes file
    ofstream out;
    out.open(outfile.c_str());
    out << total_num << " 0" << endl;
    cout << "Writing " << total_num << " electrodes to '" << outfile << "'" << endl;
    for (map_id_vector_t::iterator it = random_nodes.begin(); it != random_nodes.end(); ++it)
    {
        mat = it->first;
        for (vector_t::iterator vit = it->second->begin(); vit != it->second->end(); ++vit)
        {
            unsigned int node_id = *vit;
            out << node_id << " " << phi[node_id] << endl;
        }
    }
    out.close();

    return 0;
}
