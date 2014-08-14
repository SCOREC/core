#ifndef APF_CONVERT_H
#define APF_CONVERT_H

namespace apf {

class Mesh;
class Mesh2;
class ModelEntity;

void convert(Mesh *in, Mesh2 *out);

void construct(Mesh2* m, int* conn, int nelem, int etype);

void destruct(Mesh2* m, int*& conn, int& nelem, int &etype);

}

#endif
