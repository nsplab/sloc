/*
 * Load mesh with vcglib
 */

#include <iostream>
#include <string>
#include <sstream>
#include <boost/filesystem.hpp>
#include "color_utils.h"

// ----------------------------------------------------------------------------
// inspired by ~/dev/meshlab/src/common/meshmodel.h

#include <wrap/io_trimesh/import_stl.h>

#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/edge/base.h>
#include <vcg/simplex/face/base.h>

// OCF - optional component (fast)
#include <vcg/simplex/vertex/component_ocf.h>
#include <vcg/simplex/face/component_ocf.h>

#include <vcg/complex/used_types.h>
#include <vcg/complex/complex.h>
//#include <vcg/complex/allocate.h>

// forward declarations needed for creating the used types
class CVertexO;
class CEdgeO;
class CFaceO;

// declaration of the semantic of the used types
class CUsedTypesO : public vcg::UsedTypes<vcg::Use<CVertexO>::AsVertexType,
                                          vcg::Use<CEdgeO>::AsEdgeType,
                                          vcg::Use<CFaceO>::AsFaceType> {};


// the main vertex class
class CVertexO : public vcg::Vertex<
    CUsedTypesO,
    vcg::vertex::Coord3f,       // 12b
    vcg::vertex::BitFlags,      //  4b
    vcg::vertex::Normal3f       // 12b
> {};

// the main edge class
class CEdgeO : public vcg::Edge<
    CUsedTypesO,
    //vcg::edge::BitFlags, // 4b
    vcg::edge::EVAdj,
    vcg::edge::EEAdj
> {};

// the main face class
class CFaceO : public vcg::Face<
    CUsedTypesO,
    vcg::face::VertexRef, // 12b
    //vcg::face::BitFlags,  //  4b
    vcg::face::Normal3f,  // 12b
    //vcg::face::QualityfOcf,
    vcg::face::FFAdjOcf,
    vcg::face::VFAdjOcf
> {};

class CMeshO : public vcg::tri::TriMesh<
    vcg::vertex::vector_ocf<CVertexO>,
    vcg::face::vector_ocf<CFaceO>
> {};

// ----------------------------------------------------------------------------
namespace fs = boost::filesystem;
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cout << "Usage: "
             << ANSI_RED << argv[0] << " MESHFILE"
             << ANSI_RESET << endl;
        return 0;
    }

    fs::path path(argv[1]);
    if (!fs::exists(path))
    {
        cerr << "File "
             << ANSI_RED << path.filename()
             << ANSI_RESET << " does not exist!" << endl;
        return 1;
    }

    //MeshModel m;
    CMeshO cm;
    int loadmask;
    ostringstream errorMessage;
    int result = vcg::tri::io::ImporterSTL<CMeshO>::Open(cm, path.filename().c_str(), loadmask);
    if (result != 0) // all the importers return 0 on success
    {
        if (vcg::tri::io::ImporterSTL<CMeshO>::ErrorCritical(result))
        {
            errorMessage << "Error encountered while loading file:\n\""
                         << ANSI_RED << path.filename() << ANSI_RESET
                         << "\"\n\nError details: "
                         << tri::io::ImporterSTL<CMeshO>::ErrorMsg(result);
            cerr << errorMessage.str() << endl;
        }
    }

    return 0;
}
