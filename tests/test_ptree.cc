/*
 * Test loading of property tree from json file.
 */

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <string>
#include <vector>
#include <cassert>

struct Model
{
    std::string name;
    std::string surface_mesh;
    std::string material_data;
    std::vector<std::pair<std::string,std::string> > electrodes;
};

void read_models(std::string filename, std::vector<Model>& models)
{
    using boost::property_tree::ptree;

    ptree root;
    boost::property_tree::read_json(filename, root);

    ptree& models_pt = root.get_child("models");
    for (ptree::iterator m = models_pt.begin(); m != models_pt.end(); ++m)
    {
        Model model;

        // get the model name
        model.name = m->first;

        // value is another tree...unpack it
        ptree& model_pt = m->second;
        model.surface_mesh = model_pt.get<std::string>("mesh");
        model.material_data = model_pt.get<std::string>("material");

        // electrodes
        ptree& electrodes_pt = model_pt.get_child("electrodes");
        for (ptree::iterator e = electrodes_pt.begin(); e != electrodes_pt.end(); ++e)
        {
            model.electrodes.push_back(std::pair<std::string,std::string>(e->first, e->second.data()));
        }

        models.push_back(model);
    }

    // look for a key that's not there
    ptree::const_assoc_iterator it = root.find("blah");
    assert(it == root.not_found());

}

int main(int argc, char *argv[])
{
    if (argc > 1)
    {
        std::vector<Model> models;
        read_models(argv[1], models);
        for (std::vector<Model>::const_iterator it = models.begin(); it != models.end(); ++it)
        {
            std::cout << "model "
                      << "name='" << it->name << "' "
                      << "mesh='" << it->surface_mesh << "' "
                      << "material='" << it->material_data << "' "
                      << "electrodes=[";
            for (std::vector<std::pair<std::string,std::string> >::const_iterator j = it->electrodes.begin();
                    j != it->electrodes.end(); ++j)
            {
                std::cout << "('" << j->first << "','" << j->second << "'),";
            }
            std::cout << "]" << std::endl;
        }
    }
    else
    {
        std::cout << "Usage: " << argv[0] << " filename.json" << std::endl;
    }

    return 0;
}

