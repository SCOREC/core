/*
 * this test verifies that the convert process properly clones the underlying
 * fields, numberings, and tags
 */
#include <PCU.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <gmi_null.h>
#include <pcu_util.h>
#include <cassert>
apf::Mesh2* createEmptyMesh()
{
  gmi_model* mdl = gmi_load(".null");
  return apf::makeEmptyMdsMesh(mdl, 3, false);
}
apf::Mesh2* createMesh()
{
  apf::Mesh2* m = createEmptyMesh();
  apf::MeshEntity* verts[2];
  verts[0] =
      m->createVertex(NULL, apf::Vector3(0, 0, 0), apf::Vector3(0, 0, 0));
  verts[1] =
      m->createVertex(NULL, apf::Vector3(1, 0, 0), apf::Vector3(1, 0, 0));
  m->createEntity(apf::Mesh::Type::EDGE, NULL, verts);
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
  PCU_Comm_Init();
  gmi_register_null();
  apf::Mesh* m1 = createMesh();
  apf::Mesh2* m2 = createEmptyMesh();
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
      m1, "noField", apf::getShape(f), apf::countComponents(f));
  // create an integer tag
  apf::MeshTag* intTag = m1->createIntTag("intTag", 1);
  // loop over all vertices in mesh and set values
  int count = 1;
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
    m1->setIntTag(vert, intTag, &count);
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
  apf::Field* f2 = m2->findField("field1");
  apf::Field* uf2 = m2->findField("ufield1");

  apf::Numbering* numWithField2 = m2->findNumbering(apf::getName(numWithField));
  apf::Numbering* numNoField2 = m2->findNumbering(apf::getName(numNoField));
  apf::GlobalNumbering* globalNumWithField2 =
      m2->findGlobalNumbering(apf::getName(globalNumWithField));
  apf::GlobalNumbering* globalNumNoField2 =
      m2->findGlobalNumbering(apf::getName(globalNumNoField));
  // all of these numberings should exist
  assert(numWithField2 && numNoField2 && globalNumWithField2 &&
         globalNumNoField2);
  // make sure the fields of the numberings match up properly as they should be
  // the new field copied into the new mesh
  assert(getField(numWithField2) == f2);
  assert(getField(globalNumWithField2) == f2);
  // update the user field to reference the field in mesh 2
  apf::updateUserField(uf2, new twox(f2));
  // find the copied tag data
  apf::MeshTag* intTag2 = m2->findTag("intTag");
  assert(intTag2);
  count = 1;
  it = m2->begin(0);
  while (apf::MeshEntity* vert = m2->iterate(it)) {
    double val = apf::getScalar(f2, vert, 0);
    double uval = apf::getScalar(uf2, vert, 0);
    assert(std::abs(val - count) < 1E-15);
    assert(std::abs(uval - 2 * double(count)) < 1E-15);
    apf::setScalar(f2, vert, 0, 18.0);
    // make sure that the function updated properly
    uval = apf::getScalar(uf2, vert, 0);
    assert(std::abs(uval - 36.0) < 1E-15);
    // check to make sure the numberings have the correct values
    assert(getNumber(numWithField2, vert, 0, 0) == count);
    assert(getNumber(numNoField2, vert, 0, 0) == count);
    assert(getNumber(globalNumWithField2, vert, 0, 0) == count);
    assert(getNumber(globalNumNoField2, vert, 0, 0) == count);
    // check that the correct tag data was recovered
    int* data = new int[1];
    m2->getIntTag(vert, intTag2, data);
    assert(*data == count);
    delete[] data;
    ++count;
  }
  m2->end(it);
  delete func;
  m1->destroyNative();
  m2->destroyNative();
  apf::destroyMesh(m1);
  apf::destroyMesh(m2);
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}
