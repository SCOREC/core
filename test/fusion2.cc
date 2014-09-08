#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <parma.h>
#include <apfZoltan.h>

static apf::Mesh2* makeEmptyMesh()
{
  /* dont use null model in production */
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(g, 2, false);
  /* both planes create the field at the beginning */
  apf::createPackedField(m, "fusion", 6);
  return m;
}

static void addOneTri(apf::Mesh2* m)
{
  apf::Vector3 points[3] =
  {apf::Vector3(0,0,0),
   apf::Vector3(1,0,0),
   apf::Vector3(0,1,0)};
  apf::ModelEntity* tmp = m->findModelEntity(2,0);
  /* dont use this in production */
  apf::buildOneElement(m, tmp, apf::Mesh::TRIANGLE,points);
  apf::deriveMdsModel(m);
}

static void setValues(apf::Mesh2* m)
{
  apf::Field* f = m->findField("fusion");
  apf::MeshEntity* e = apf::getMdsEntity(m, 2, 0);
  apf::MeshEntity* v[3];
  m->getDownward(e, 0, v);
  double values[6] = {1.,2.,3.,4.,5.,6.};
  for (int i = 0; i < 3; ++i)
    apf::setComponents(f, v[i], 0, values);
  /* dont use this either */
}

static void ensureNewNumbering(apf::Mesh2* m)
{
  apf::FieldShape* s = m->getShape();
  apf::Numbering* local = m->findNumbering(s->getName());
  if (local) {
    fprintf(stderr,"removing old numbering\n");
    apf::destroyNumbering(local);
  }
}

static void storeInArray(apf::Mesh2* m)
{
  ensureNewNumbering(m);
  /* remember findField is a linear search,
     in production keep the pointer around */
  apf::freeze(m->findField("fusion"));
}

static void backToTags(apf::Mesh2* m)
{
  apf::unfreeze(m->findField("fusion"));
}

static void connectPlanes(apf::Mesh2* m)
{
  /* dont use any of this function in production. */
  int peer = 1 - PCU_Comm_Self();
  for (int i = 0; i < 3; ++i) {
    apf::MeshEntity* v = apf::getMdsEntity(m, 0, i);
    m->addRemote(v, peer, v); //pointers aren't usually equal
  }
  apf::stitchMesh(m);
  m->acceptChanges();
}

static void checkValues(apf::Mesh2* m)
{
  apf::Field* f = m->findField("fusion");
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it)))
    assert( apf::hasEntity(f, v) );
  m->end(it);
}

struct GroupCode : public Parma_GroupCode
{
  apf::Mesh2* mesh;
  void run(int group)
  {
    mesh = ::makeEmptyMesh();
    if (group == 0) {
      addOneTri(mesh);
      setValues(mesh);
      checkValues(mesh);
      storeInArray(mesh);
/* this is where m3dc1 would run the 2D solver */
      backToTags(mesh);
      checkValues(mesh);
    }
  }
};

static void copyPlane(apf::Mesh2* m)
{
  /* this mimics the copy entity construction */
  if (PCU_Comm_Self() == 1)
    addOneTri(m);
  /* this connects the two one-triangle planes together */
  connectPlanes(m);
}

/* send field values from rank 0 to rank 1.
   we are just lucky that lower ranks win ownership
   in this 2-rank test case.
   use a Sharing object in production*/
static void syncValues(apf::Mesh2* m)
{
  apf::synchronize(m->findField("fusion"));
}

static void globalCode(apf::Mesh2* m)
{
  copyPlane(m);
  syncValues(m);
  checkValues(m);
  storeInArray(m);
  checkValues(m);
  /* run the 3D solve... */
  backToTags(m);
  checkValues(m);
  m->destroyNative();
  apf::destroyMesh(m);
}

int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  assert(PCU_Comm_Peers() == 2);
  gmi_register_null();
  GroupCode code;
  int const groupSize = 1;
  apf::Unmodulo outMap(PCU_Comm_Self(), groupSize);
  Parma_SplitPartition(NULL, groupSize, code);
  /* update mesh for leaving groups */
  apf::remapPartition(code.mesh, outMap);
  globalCode(code.mesh);
  PCU_Comm_Free();
  MPI_Finalize();
}
