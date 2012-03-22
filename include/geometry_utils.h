/* 
 * Collection of useful geometry functions
 */
#ifndef SLOC_GEOMETRY_UTILS_H
#define SLOC_GEOMETRY_UTILS_H

#include "mesh.h"

namespace sloc
{
// ----------------------------------------------------------------------------

void centroid(double A[3], double B[3], double C[3], double D[3]);
void midpoint(double A[3], double B[3], double M[3]);

/* Convert triangles into quadrilateral elements. */
void tri2quad(long t[3], long u[4], long q[3*4]);
void tri2quad(Mesh& tmesh, Mesh& qmesh);

// ----------------------------------------------------------------------------
}
#endif
