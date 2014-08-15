#ifndef APF_CONVERT_H
#define APF_CONVERT_H

#include <map>

namespace apf {

class Mesh;
class Mesh2;
class ModelEntity;
class MeshEntity;

void convert(Mesh *in, Mesh2 *out);

typedef std::map<int, MeshEntity*> GlobalToVert;

void construct(Mesh2* m, int* conn, int nelem, int etype,
    GlobalToVert& globalToVert);

/** \brief Assign coordinates to the mesh
  * Each peer provides a set of the coordinates. The coords most be ordered
  * according to the global ids of the vertices. Peer 0 provides the coords
  * for vertices 0 to m-1, peer to for m to n-1, ...
  */
void setCoords(Mesh2* m, const double* coords, int nverts,
    const GlobalToVert& globalToVert);

void destruct(Mesh2* m, int*& conn, int& nelem, int &etype);

void extractCoords(Mesh2* m, double*& coords, int& nverts);

}

#endif
