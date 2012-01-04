/* Convert triangles into quadrilateral elements.
 *
 */
#ifndef TRI2QUAD_H
#define TRI2QUAD_H

#include "mesh.h"

namespace sloc
{
// ----------------------------------------------------------------------------

void centroid(double A[3], double B[3], double C[3], double D[3]);
void midpoint(double A[3], double B[3], double M[3]);
void tri2quad(long t[3], long u[4], long q[3*4]);
void tri2quad(Mesh& tmesh, Mesh& qmesh);

// ----------------------------------------------------------------------------
}
#endif
