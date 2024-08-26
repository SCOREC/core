/*
 * this test verifies that the convert process properly clones the underlying
 * fields, numberings, and tags
 */
#include <lionPrint.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <gmi_null.h>
#include <pcu_util.h>
#include <cassert>

apf::Mesh2* createEmptyMesh(pcu::PCU *PCUObj)
{
  gmi_model* mdl = gmi_load(".null");
  return apf::makeEmptyMdsMesh(mdl, 1, false, PCUObj);
}
apf::Mesh2* createMesh(pcu::PCU *PCUObj)
{
  apf::Mesh2* m = createEmptyMesh(PCUObj);
  apf::Vector3 pts[2] = {apf::Vector3(0,0,0), apf::Vector3(1,0,0)};
  apf::MeshEntity* verts[2];
  for( int i=0; i<2; i++)
    verts[i] = m->createVertex(m->findModelEntity(0,i), pts[i], pts[i]);
  apf::buildElement(m, m->findModelEntity(1,2), apf::Mesh::EDGE,verts);
  m->acceptChanges();
  m->verify();
  return m;
}
class twox : public apf::Function {
  private:
  apf::Field* x;

  public:
  virtual void eval(apf::MeshEntity* e, double* result)
  {
    result[0] = 2 * apf::getScalar(x, e, 0);
  }
  twox(apf::Field* x) : x(x) {}
};
int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_null();

  // create meshes and write data to one of them
  apf::Mesh* m1 = createMesh(&PCUObj);
  apf::Mesh2* m2 = createEmptyMesh(&PCUObj);
  // create field on m1
  apf::Field* f = apf::createLagrangeField(m1, "field1", apf::SCALAR, 1);
  apf::Function* func = new twox(f);
  apf::Field* uf =
      apf::createUserField(m1, "ufield1", apf::SCALAR, apf::getShape(f), func);
  // create numbering and global numbering with and without fields to make sure
  // they transfer properly
  apf::Numbering* numWithField = apf::createNumbering(f);
  apf::Numbering* numNoField = apf::createNumbering(
      m1, "noField", apf::getShape(f), apf::countComponents(f));
  apf::GlobalNumbering* globalNumWithField = apf::createGlobalNumbering(f);

  apf::GlobalNumbering* globalNumNoField = apf::createGlobalNumbering(
      m1, "noField_global", apf::getShape(f), apf::countComponents(f));
  // create an integer tag
  apf::MeshTag* intTag = m1->createIntTag("intTag", 1);
  // loop over all vertices in mesh and set values
  long count = 1;
  apf::MeshIterator* it = m1->begin(0);
  while (apf::MeshEntity* vert = m1->iterate(it)) {
    apf::setScalar(f, vert, 0, count);
    // set the numberings
    apf::number(numWithField, vert, 0, 0, count);
    apf::number(numNoField, vert, 0, 0, count);
    apf::number(globalNumWithField, vert, 0, count);
    apf::number(globalNumNoField, vert, 0, count);
    // set the tag
    // int tagData[1] = {count};
    int v = static_cast<int>(count);
    m1->setIntTag(vert, intTag, &v);
    ++count;
  }
  m1->end(it);
  // verify the user field works on mesh 1
  it = m1->begin(0);
  while (apf::MeshEntity* vert = m1->iterate(it)) {
    double val = apf::getScalar(f, vert, 0);
    double uval = apf::getScalar(uf, vert, 0);
    if (!(std::abs(uval - 2 * double(val)) < 1E-15)) return 1;
  }
  m1->end(it);


  // copy m1 to m2
  apf::convert(m1, m2);
  m2->verify();
  apf::Field* f2 = m2->findField("field1");
  apf::Field* uf2 = m2->findField("ufield1");

  // all of these numberings should exist
  assert(m2->findNumbering(apf::getName(numWithField)));
  assert(m2->findNumbering(apf::getName(numNoField)));
  assert(m2->findGlobalNumbering(apf::getName(globalNumWithField)));
  assert(m2->findGlobalNumbering(apf::getName(globalNumNoField)));
  // make sure the fields of the numberings match up properly as they should be
  // the new field copied into the new mesh
  assert(getField(m2->findNumbering(apf::getName(numWithField))) == f2);
  assert(getField(m2->findGlobalNumbering(apf::getName(globalNumWithField))) == f2);
  // update the user field to reference the field in mesh 2
  apf::Function* func2 = new twox(f2);
  apf::updateUserField(uf2, func2);
  // find the copied tag data
  apf::MeshTag* intTag2 = m2->findTag("intTag");
  assert(intTag2);
  count = 1;
  it = m2->begin(0);
  while (apf::MeshEntity* vert = m2->iterate(it)) {
    assert(std::abs(apf::getScalar(f2, vert, 0) - count) < 1E-15);
    assert(std::abs(apf::getScalar(uf2, vert, 0) - 2 * double(count)) < 1E-15);
    apf::setScalar(f2, vert, 0, 18.0);
    // make sure that the function updated properly
    assert(std::abs(apf::getScalar(uf2, vert, 0) - 36.0) < 1E-15);
    // check to make sure the numberings have the correct values
    assert(getNumber(m2->findNumbering(apf::getName(numWithField)), vert, 0, 0) == count);
    assert(getNumber(m2->findNumbering(apf::getName(numNoField)), vert, 0, 0) == count);
    assert(getNumber(m2->findGlobalNumbering(apf::getName(globalNumWithField)), vert, 0, 0) == count);
    assert(getNumber(m2->findGlobalNumbering(apf::getName(globalNumNoField)), vert, 0, 0) == count);
    // check that the correct tag data was recovered
    int* data = new int[1];
    m2->getIntTag(vert, intTag2, data);
    assert(*data == count);
    delete[] data;
    ++count;
  }
  m2->end(it);

  // check that not transfering Fields/Numberings/Tags also works
  apf::Mesh2* m3 = createEmptyMesh(&PCUObj);
  apf::convert(m1, m3, NULL, NULL, false);
  m3->verify();

  assert(!m3->findNumbering(apf::getName(numWithField)));
  assert(!m3->findNumbering(apf::getName(numNoField)));
  assert(!m3->findGlobalNumbering(apf::getName(globalNumWithField)));
  assert(!m3->findGlobalNumbering(apf::getName(globalNumNoField)));
  assert(!m3->findNumbering(apf::getName(numWithField)));
  assert(!m3->findGlobalNumbering(apf::getName(globalNumWithField)));
  assert(!m3->findTag("intTag"));



  // cleanup
  delete func;
  delete func2;
  m1->destroyNative();
  m2->destroyNative();
  m3->destroyNative();
  apf::destroyMesh(m1);
  apf::destroyMesh(m2);
  apf::destroyMesh(m3);
  }
  MPI_Finalize();
  return 0;
}
