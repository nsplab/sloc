#ifndef IO_UCD_H
#define IO_UCD_H

#include <string>
#include <vector>
#include <valarray>
#include <fstream>

namespace sloc
{
// ----------------------------------------------------------------------------

typedef enum UCD_Cell_Type
{
    UCD_CELL_NONE=0,
    UCD_CELL_TRI,
    UCD_CELL_QUAD
} UCD_Cell_Type;

UCD_Cell_Type string2celltype(const std::string& cell_type);
std::string celltype2string(UCD_Cell_Type cell_type);

class UCD_Cell
{
public:
    UCD_Cell();
    ~UCD_Cell();
    int num_vertices() const;
public:
    long cell_id;
    long mat_id;
    UCD_Cell_Type cell_type;
    std::vector<long> cell_verts;
};

// ----------------------------------------------------------------------------

class UCD_File
{
public:
    UCD_File();
    ~UCD_File();

    void clear();

    void read(const char *filename);
    void write(const char *filename);

    void _read_header(std::ifstream& is);
    void _read_points(std::ifstream& is);
    void _read_cells(std::ifstream& is);

    void _write_header(std::ofstream& os);
    void _write_points(std::ofstream& os);
    void _write_cells(std::ofstream& os);

public:
    long num_nodes;
    long num_cells;
    long num_ndata;
    long num_cdata;
    long num_mdata;
    std::vector<std::string> _comments;
    std::valarray<double> _nodes;
    std::valarray<long> _node_ids;
    std::vector<UCD_Cell*> _cells;
};

// ----------------------------------------------------------------------------
}

#endif
