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

void destruct(Mesh2* m, int*& conn, int& nelem, int &etype);

}

#endif
