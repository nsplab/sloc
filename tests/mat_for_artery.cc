/* artery.cc
 *
 * Here, we assign a different material id to the quad elements
 * in Artery2.quad.ucd
 */

#include <iostream>
#include "io_ucd.h"

int main(void)
{
    std::string infile = "tmp/Artery2.quad.ucd";
    std::string outfile = "tmp/Artery2_1.quad.ucd";

    sloc::UCD_File ucd;
    ucd.read(infile.c_str());
    std::cout << "Read " << infile << std::endl;
    
    // Since Bone2.quad.ucd uses material id 1,
    // we will use material id 2 for elements in
    // Artery2.quad.ucd
    long mat_id = 2;
    std::vector<sloc::UCD_Cell*>::iterator it;
    for (it = ucd._cells.begin(); it != ucd._cells.end(); ++it)
    {
        sloc::UCD_Cell *cell = *it;
        cell->mat_id = mat_id;
    }

    ucd.write(outfile.c_str());
    std::cout << "Wrote " << outfile << std::endl;

    return 0;
}
