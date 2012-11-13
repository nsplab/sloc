/*
 * The following two links are a good reference for the UCD format:
 *   http://people.sc.fsu.edu/~jburkardt/html/ucd_format.html
 *   http://people.sc.fsu.edu/~jburkardt/data/ucd/ucd.html
 */

#include "sloc/io_ucd.h"
#include <iostream>
#include <cassert>

using namespace std;
using namespace sloc;

// ----------------------------------------------------------------------------

UCD_Cell::UCD_Cell()
    : cell_id(-1),
      mat_id(-1),
      cell_type(UCD_CELL_NONE)
{
}

UCD_Cell::~UCD_Cell()
{
}

int UCD_Cell::num_vertices() const
{
    switch (cell_type)
    {
    case UCD_CELL_TRI:
        return 3;
    case UCD_CELL_QUAD:
        return 4;
    default:
        return 0;
    }
}

UCD_Cell_Type sloc::string2celltype(const string& cell_type)
{
    if (cell_type == "tri")
        return UCD_CELL_TRI;
    if (cell_type == "quad")
        return UCD_CELL_QUAD;
    return UCD_CELL_NONE;
}

string sloc::celltype2string(UCD_Cell_Type cell_type)
{
    switch (cell_type)
    {
    case UCD_CELL_TRI:
        return "tri";
    case UCD_CELL_QUAD:
        return "quad";
    default:
        return "none";
    }
}

// ----------------------------------------------------------------------------

UCD_File::UCD_File()
    : num_nodes(0),
      num_cells(0),
      num_ndata(0),
      num_cdata(0),
      num_mdata(0)
{
}

UCD_File::~UCD_File()
{
    clear();
}

void UCD_File::clear()
{
    _comments.clear();
    num_nodes = 0;
    num_cells = 0;
    num_ndata = 0;
    num_cdata = 0;
    num_mdata = 0;

    _nodes.resize(0);
    _node_ids.resize(0);

    for (unsigned int i = 0; i < _cells.size(); i++)
        delete _cells[i];
    _cells.clear();
}

// ----------------------------------------------------------------------------

void UCD_File::read(const char *filename)
{
    ifstream is;
    is.open(filename);
    clear();
    _read_header(is);
    _read_points(is);
    _read_cells(is);
    is.close();
}

void UCD_File::write(const char *filename)
{
    ofstream os;
    os.open(filename);
    _write_header(os);
    _write_points(os);
    _write_cells(os);
    os.close();
}

// ----------------------------------------------------------------------------

void UCD_File::_read_header(ifstream& is)
{
    string line;
    while (is.peek() == '#')
    {
        getline(is, line);
        _comments.push_back(line);
    }

    is >> num_nodes;
    is >> num_cells;
    is >> num_ndata;
    is >> num_cdata;
    is >> num_mdata;

    _nodes.resize(num_nodes * 3);
    _node_ids.resize(num_nodes);
    _cells.reserve(num_cells);
}

void UCD_File::_read_points(ifstream& is)
{
    assert(num_nodes > 0);

    long node_id;
    double x, y, z;

    for (int n = 0; n < num_nodes; n++)
    {
        is >> node_id >> x >> y >> z;
        _node_ids[n] = node_id;
        _nodes[3*n+0] = x;
        _nodes[3*n+1] = y;
        _nodes[3*n+2] = z;
    }
}

void UCD_File::_read_cells(ifstream& is)
{
    assert(num_cells > 0);

    string s;
    int nno;
    long id;

    for (long e = 0; e < num_cells; e++)
    {
        UCD_Cell *cell = new UCD_Cell;

        is >> (cell->cell_id);

        is >> (cell->mat_id);

        is >> s;
        cell->cell_type = string2celltype(s);

        nno = cell->num_vertices();
        for (int i = 0; i < nno; i++)
        {
            is >> id;
            cell->cell_verts.push_back(id);
        }

        _cells.push_back(cell);
    }

}

// ----------------------------------------------------------------------------

void UCD_File::_write_header(ofstream& os)
{
    vector<string>::iterator it;
    for (it = _comments.begin(); it != _comments.end(); ++it)
        os << *it << endl;

    os << num_nodes << " "
       << num_cells << " "
       << num_ndata << " "
       << num_cdata << " "
       << num_mdata << endl;
}

void UCD_File::_write_points(ofstream& os)
{
    for (int n = 0; n < num_nodes; n++)
    {
        os << _node_ids[n] << " "
           << _nodes[3*n+0] << " "
           << _nodes[3*n+1] << " "
           << _nodes[3*n+2] << endl;
    }
}

void UCD_File::_write_cells(ofstream& os)
{
    for (long e = 0; e < num_cells; e++)
    {
        UCD_Cell *cell = _cells[e];

        os << (cell->cell_id) << " ";
        os << (cell->mat_id) << " ";

        string cell_type = celltype2string(cell->cell_type);
        os << cell_type;

        int nno = cell->num_vertices();
        for (int i = 0; i < nno; i++)
            os << " " << cell->cell_verts[i];

        os << endl;
    }
}

// ----------------------------------------------------------------------------
// EOF
