#include "crv.h"
#include "crvQuality.h"
#include "crvDBG.h"
#include <gmi.h>


namespace crv_dbg
{

void printTetNumber(apf::Mesh2* m, apf::MeshEntity* e, const char* numberingName)
{
  apf::Numbering* n = m->findNumbering(numberingName);
  if (!n) return;
  printf("TET:: %d\n", apf::getNumber(n, e, 0, 0));
}

void printEdgeCavityInvalidities(apf::Mesh2* m, apf::MeshEntity* e, apf::Numbering* n)
{
  if (n)
    printf("entity %d, ", apf::getNumber(n, e, 0, 0));
  else
    printf("entity, ");
  printf("on model %d\n", m->getModelTag(m->toModel(e)));
  apf::Adjacent upTets;
  m->getAdjacent(e, 3, upTets);
  for (std::size_t i = 0; i < upTets.getSize(); i++) {
    std::vector<int> ai = crv::getAllInvalidities(m, upTets[i]);
    for (std::size_t j = 0; j < ai.size(); j++) {
      printf("%d ", ai[j]);
    }
    printf("\n");
  }
}

void visualizeCavityMesh(apf::Mesh2* m, apf::MeshEntity* ent,
    const char* prefix, apf::Numbering* n, int resolution)
{
  int num = 0;
  if (n)
    num = apf::getNumber(n, ent, 0, 0);

  int dim = m->getDimension();
  apf::Adjacent e;
  m->getAdjacent(ent, dim, e);


  gmi_model* g = gmi_load(".null");
  apf::Mesh2* outMesh = apf::makeEmptyMdsMesh(g, dim, false);

  apf::MeshEntity* newEnt[99];

  for (std::size_t ii = 0; ii < e.getSize(); ii++) {
    // Verts
    apf::MeshEntity* vs[12];
    apf::MeshEntity* newVs[12];
    int nv = m->getDownward(e[ii], 0, vs);
    for(int i = 0; i < nv; ++i)
    {
      apf::Vector3 p;
      apf::Vector3 param(0., 0., 0.);
      m->getPoint(vs[i], 0, p);
      newVs[i] = outMesh->createVertex(0, p, param);
    }

    // Edges
    apf::MeshEntity* es[12];
    apf::MeshEntity* newEs[12];
    int ne = m->getDownward(e[ii], 1, es);
    for(int i = 0; i < ne; ++i)
    {
      apf::MeshEntity* evs[2];
      apf::MeshEntity* new_evs[2];
      m->getDownward(es[i], 0, evs);
      for (int j = 0; j < 2; j++) {
        new_evs[j] = newVs[apf::findIn(vs, nv, evs[j])];
      }

      newEs[i] = outMesh->createEntity(apf::Mesh::EDGE, 0, new_evs);
    }

    // Faces
    apf::MeshEntity* fs[12];
    apf::MeshEntity* newFs[12];
    int nf = m->getDownward(e[ii], 2, fs);
    for(int i = 0; i < nf; ++i)
    {
      apf::MeshEntity* fes[3];
      apf::MeshEntity* new_fes[3];
      m->getDownward(fs[i], 1, fes);
      for (int j = 0; j < 3; j++) {
        new_fes[j] = newEs[apf::findIn(es, ne, fes[j])];
      }
      newFs[i] = outMesh->createEntity(apf::Mesh::TRIANGLE, 0, new_fes);
    }

    // Regions
    apf::MeshEntity* tet = 0;
    if (dim == 3) {
      tet = outMesh->createEntity(apf::Mesh::TET, 0, newFs);
    }

    newEnt[ii] = dim == 2 ? newFs[0] : tet;

    PCU_ALWAYS_ASSERT(m->getType(e[ii]) == outMesh->getType(newEnt[ii]));
    outMesh->acceptChanges();
  }
  apf::changeMeshShape(outMesh, crv::getBezier(3), true);
  outMesh->acceptChanges();

  for (std::size_t ii = 1; ii < e.getSize(); ii++) {
    for (int d = 1; d <= dim; d++)
    {
      if (!m->getShape()->hasNodesIn(d)) continue;
      apf::MeshEntity* eds[12];
      int counter = m->getDownward(e[ii], d, eds);
      apf::MeshEntity* new_eds[12];
      outMesh->getDownward(newEnt[ii], d, new_eds);
      int non = outMesh->getShape()->countNodesOn(apf::Mesh::simplexTypes[d]);
      for(int n = 0; n < counter; ++n) {
        for(int i = 0; i < non; ++i) {
        apf::Vector3 p;
        m->getPoint(eds[n], i, p);
        outMesh->setPoint(new_eds[n], i, p);
        }
      }
    }
    outMesh->acceptChanges();
  }

  std::stringstream ss;
  ss << prefix<< num;
  crv::writeCurvedVtuFiles(outMesh, apf::Mesh::TET, resolution, ss.str().c_str());
  crv::writeCurvedVtuFiles(outMesh, apf::Mesh::TRIANGLE, resolution, ss.str().c_str());
  crv::writeCurvedWireFrame(outMesh, resolution, ss.str().c_str());

  outMesh->destroyNative();
  apf::destroyMesh(outMesh);
}


void visualizeIndividualCavityEntities(apf::Mesh2* m, apf::MeshEntity* ent,
    const char* prefix, apf::Numbering* n, int resolution)
{
  int num = 0;
  if (n)
    num = apf::getNumber(n, ent, 0, 0);

  int dim = m->getDimension();
  apf::Adjacent e;
  m->getAdjacent(ent, dim, e);


  gmi_model* g = gmi_load(".null");

  int nat = e.getSize();
  apf::Mesh2* outMesh[nat];
  for (int i = 0; i < nat; i++) {
    outMesh[i] = apf::makeEmptyMdsMesh(g, dim, false);
  }

  apf::MeshEntity* newEnt[nat];

  for (int ii = 0; ii < nat; ii++) {
    // Verts
    apf::MeshEntity* vs[12];
    apf::MeshEntity* newVs[12];
    int nv = m->getDownward(e[ii], 0, vs);
    for(int i = 0; i < nv; ++i)
    {
      apf::Vector3 p;
      apf::Vector3 param(0., 0., 0.);
      m->getPoint(vs[i], 0, p);
      newVs[i] = outMesh[ii]->createVertex(0, p, param);
    }

    // Edges
    apf::MeshEntity* es[12];
    apf::MeshEntity* newEs[12];
    int ne = m->getDownward(e[ii], 1, es);
    for(int i = 0; i < ne; ++i)
    {
      apf::MeshEntity* evs[2];
      apf::MeshEntity* new_evs[2];
      m->getDownward(es[i], 0, evs);
      for (int j = 0; j < 2; j++) {
        new_evs[j] = newVs[apf::findIn(vs, nv, evs[j])];
      }

      newEs[i] = outMesh[ii]->createEntity(apf::Mesh::EDGE, 0, new_evs);
    }
   // Faces
    apf::MeshEntity* fs[12];
    apf::MeshEntity* newFs[12];
    int nf = m->getDownward(e[ii], 2, fs);
    for(int i = 0; i < nf; ++i)
    {
      apf::MeshEntity* fes[3];
      apf::MeshEntity* new_fes[3];
      m->getDownward(fs[i], 1, fes);
      for (int j = 0; j < 3; j++) {
        new_fes[j] = newEs[apf::findIn(es, ne, fes[j])];
      }
      newFs[i] = outMesh[ii]->createEntity(apf::Mesh::TRIANGLE, 0, new_fes);
    }

    // Regions
    apf::MeshEntity* tet = 0;
    if (dim == 3) {
      tet = outMesh[ii]->createEntity(apf::Mesh::TET, 0, newFs);
    }

    newEnt[ii] = dim == 2 ? newFs[0] : tet;

    PCU_ALWAYS_ASSERT(m->getType(e[ii]) == outMesh[ii]->getType(newEnt[ii]));
    outMesh[ii]->acceptChanges();

    apf::changeMeshShape(outMesh[ii], crv::getBezier(m->getShape()->getOrder()), true);
    outMesh[ii]->acceptChanges();


    for (int d = 1; d <= dim; d++)
    {
      if (!m->getShape()->hasNodesIn(d)) continue;
      apf::MeshEntity* eds[12];
      int counter = m->getDownward(e[ii], d, eds);
      apf::MeshEntity* new_eds[12];
      outMesh[ii]->getDownward(newEnt[ii], d, new_eds);
      int non = outMesh[ii]->getShape()->countNodesOn(apf::Mesh::simplexTypes[d]);
      for(int n = 0; n < counter; ++n) {
        for(int i = 0; i < non; ++i) {
        apf::Vector3 p;
        m->getPoint(eds[n], i, p);
        outMesh[ii]->setPoint(new_eds[n], i, p);
        }
      }
    }
    outMesh[ii]->acceptChanges();

    std::stringstream ss;
    ss << prefix<< num << "_TET_"<< ii;
    crv::writeCurvedVtuFiles(outMesh[ii], apf::Mesh::TET, resolution, ss.str().c_str());
    crv::writeCurvedVtuFiles(outMesh[ii], apf::Mesh::TRIANGLE, resolution, ss.str().c_str());
    crv::writeCurvedWireFrame(outMesh[ii], resolution, ss.str().c_str());

    outMesh[ii]->destroyNative();
    apf::destroyMesh(outMesh[ii]);
  }
}

void visualizeTetFaces(apf::Mesh2* m, apf::MeshEntity* e, const char* prefix, int resolution)
{
  if (m->getType(e) != apf::Mesh::TET)
    return;
  int dim = 3;

  gmi_model* g = gmi_load(".null");

  apf::Mesh2* outMesh[4];
  for (int i = 0; i < 4; i++) {
    outMesh[i] = apf::makeEmptyMdsMesh(g, 2, false);
  }

  apf::MeshEntity* face[4];
  apf::MeshEntity* newface[4];
  int nf = m->getDownward(e, 2, face);

  for (int i = 0; i < nf; i++) {

    //Verts
    apf::MeshEntity* vs[3];
    apf::MeshEntity* newVs[3];
    int nv = m->getDownward(face[i], 0, vs);
    for (int j = 0; j < nv; j++) {
      apf::Vector3 p;
      apf::Vector3 param(0., 0., 0.);
      m->getPoint(vs[j], 0, p);
      newVs[j] = outMesh[i]->createVertex(0, p, param);
    }

    //Edges
    apf::MeshEntity* es[3];
    apf::MeshEntity* newEs[3];
    int ne = m->getDownward(face[i], 1, es);
    for (int j = 0; j < ne; j++) {
      apf::MeshEntity* evs[2];
      apf::MeshEntity* newEvs[2];
      m->getDownward(es[j], 0, evs);
      for (int k = 0; k < 2; k++) {
        int kk = apf::findIn(vs,nv, evs[k]);
        newEvs[k] = newVs[kk];
      }
      newEs[j] = outMesh[i]->createEntity(apf::Mesh::EDGE, 0, newEvs);
    }

    //Faces
    apf::MeshEntity* fes[3];
    apf::MeshEntity* newFes[3];
    m->getDownward(face[i], 1, fes);
    for (int j = 0; j < 3; j++)
      newFes[j] = newEs[apf::findIn(es, ne, fes[j])];

    newface[i] = outMesh[i]->createEntity(apf::Mesh::TRIANGLE, 0, newFes);

    PCU_ALWAYS_ASSERT(m->getType(face[i]) == outMesh[i]->getType(newface[i]));

    outMesh[i]->acceptChanges();

    apf::changeMeshShape(outMesh[i], crv::getBezier(m->getShape()->getOrder()), true);
    outMesh[i]->acceptChanges();

    for (int d = 1; d < dim; d++) {
      if (!m->getShape()->hasNodesIn(d)) continue;
      apf::MeshEntity* ent[12];
      int counter = m->getDownward(face[i], d, ent);
      apf::MeshEntity* newent[12];
      outMesh[i]->getDownward(newface[i], d, newent);
      int non = outMesh[i]->getShape()->countNodesOn(apf::Mesh::simplexTypes[d]);
      for (int n = 0; n < counter; n++) {
        for (int j = 0; j < non; j++) {
          apf::Vector3 p;
          m->getPoint(ent[n], j, p);
          outMesh[i]->setPoint(newent[n], j, p);
        }
      }
    }

    outMesh[i]->acceptChanges();

    std::stringstream ss;
    ss << prefix << "_Face_"<< i;
    crv::writeCurvedVtuFiles(outMesh[i], apf::Mesh::TRIANGLE, resolution, ss.str().c_str());
    crv::writeCurvedWireFrame(outMesh[i], resolution, ss.str().c_str());
  }
}


}
