#ifndef APF_MPAS_H
#define APF_MPAS_H

namespace apf {

class Mesh2;

void loadMpasMesh(apf::Mesh2* m, const char* filename);
void writeMpasAssignments(apf::Mesh2* m, const char* ncFilename,
    const char* outPrefix);

}

#endif
