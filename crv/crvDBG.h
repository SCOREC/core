#ifndef CRVDBG_H
#define CRVDBG_H

#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfNumbering.h>

namespace crv_dbg
{
void printTetNumber(apf::Mesh2* m, apf::MeshEntity* e, const char* numberingName);
void printCavityInvalidities(apf::Mesh2* m, apf::MeshEntity* e,
    apf::Numbering* n = 0);

void visualizeCavityMesh(apf::Mesh2* m, apf::MeshEntity* ent,
    const char* prefix, apf::Numbering* n = 0, int resolution = 8);

void visualizeIndividualCavityEntities(apf::Mesh2* m, apf::MeshEntity* ent,
    const char* prefix, apf::Numbering* n = 0, int resolution = 8);

void visualizeTetFaces(apf::Mesh2* m, apf::MeshEntity* e,
    const char* prefix, int resolution = 8);
}
#endif
