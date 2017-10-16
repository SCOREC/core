/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include "apfNumbering.h"
#include "apfNumberingClass.h"
#include "apfField.h"
#include "apfMesh.h"
#include "apfShape.h"
#include "apfTagData.h"
#include <pcu_util.h>

namespace apf {

enum { FIXED = -2, FREE_BUT_NOT_NUMBERED = -1 };

template <class T>
NumberingOf<T>::NumberingOf()
{
  field = 0;
  components = 0;
}

template <class T>
int NumberingOf<T>::countComponents() const {return components;}

template <>
int NumberingOf<int>::getScalarType() {return Mesh::INT;}
template <>
int NumberingOf<long>::getScalarType() {return Mesh::LONG;}

template <class T>
void NumberingOf<T>::init(
    const char* n,
    Mesh* m,
    FieldShape* s,
    int c)
{
  components = c;
  FieldBase::init(n,m,s,new TagDataOf<T>());
}

template <class T>
void NumberingOf<T>::init(Field* f)
{
  field = f;
  std::string nm = f->getName();
  nm += "_num";
  init(nm.c_str(),
      f->getMesh(),
      f->getShape(),
      f->countComponents());
}

template <class T>
Field* NumberingOf<T>::getField() {return field;}

template <class T>
FieldDataOf<T>* NumberingOf<T>::getData()
{
  return static_cast<FieldDataOf<T>*>(FieldBase::getData());
}

template <class T>
void NumberingOf<T>::getAll(MeshEntity* e, T* dat)
{
  FieldDataOf<T>* fieldData = getData();
  if (fieldData->hasEntity(e))
    fieldData->get(e,dat);
  else
  { //default initialization to free and not numbered
    int n = countValuesOn(e);
    for (int i=0; i < n; ++i)
      dat[i] = FREE_BUT_NOT_NUMBERED;
  }
}

template <class T>
T NumberingOf<T>::get(MeshEntity* e, int node, int component)
{
  NewArray<T> data(countValuesOn(e));
  getAll(e,&(data[0]));
  return data[node*components + component];
}

template <class T>
void NumberingOf<T>::set(MeshEntity* e, int node, int component, T value)
{
  NewArray<T> data(countValuesOn(e));
  getAll(e,&(data[0]));
  data[node*components + component] = value;
  getData()->set(e,&(data[0]));
}

/* explicit instantiations */
template class NumberingOf<int>;
template class NumberingOf<long>;

Numbering* createNumbering(Field* f)
{
  Numbering* n = new Numbering();
  n->init(f);
  f->getMesh()->addNumbering(n);
  return n;
}

Numbering* createNumbering(
    Mesh* mesh,
    const char* name,
    FieldShape* shape,
    int components)
{
  Numbering* n = new Numbering();
  n->init(name,mesh,shape,components);
  mesh->addNumbering(n);
  return n;
}

void destroyNumbering(Numbering* n)
{
  n->getMesh()->removeNumbering(n);
  delete n;
}

void fix(Numbering* n, MeshEntity* e, int node, int component, bool fixed)
{
  if (fixed)
    n->set(e,node,component,FIXED);
  else
    n->set(e,node,component,FREE_BUT_NOT_NUMBERED);
}

bool isFixed(Numbering* n, MeshEntity* e, int node, int component)
{
  return FIXED == n->get(e,node,component);
}

bool isNumbered(Numbering* n, MeshEntity* e, int node, int component)
{
  return n->get(e,node,component) >= 0;
}

void number(Numbering* n, MeshEntity* e, int node, int component, int number)
{
  PCU_ALWAYS_ASSERT( ! isFixed(n,e,node,component));
  PCU_ALWAYS_ASSERT(number >= 0);
  n->set(e,node,component,number);
}

int getNumber(Numbering* n, MeshEntity* e, int node, int component)
{
  PCU_ALWAYS_ASSERT(isNumbered(n,e,node,component));
  return n->get(e,node,component);
}

Field* getField(Numbering* n)
{
  return n->getField();
}

FieldShape* getShape(Numbering* n)
{
  return n->getShape();
}

const char* getName(Numbering* n)
{
  return n->getName();
}

Mesh* getMesh(Numbering* n)
{
  return n->getMesh();
}

int countComponents(Numbering* n)
{
  return n->countComponents();
}

int getElementNumbers(Numbering* n, MeshEntity* e, NewArray<int>& numbers)
{
  return n->getData()->getElementData(e,numbers);
}

int countFixed(Numbering* n)
{
  int num_fixed = 0;
  Mesh* m = n->getMesh();
  FieldShape* fs = n->getShape();
  int num_components = n->getField()->countComponents();
  
  for(int ii = 0; ii < 3 && fs->hasNodesIn(ii); ii++)
  {
    apf::MeshIterator * it = m->begin(ii);
    while(apf::MeshEntity * en = m->iterate(it))
    {
      if(m->isOwned(en))
      {
        for(int jj = 0; jj < fs->countNodesOn(m->getType(en)); jj++)
          for(int kk = 0; kk < num_components; kk++)
            num_fixed = isFixed(n,en,jj,kk) ? num_fixed + 1 : num_fixed;
      }
    }
  }
  return num_fixed;
}

void synchronize(Numbering * n, Sharing* shr)
{
  synchronizeFieldData<int>(n->getData(), shr);
}

struct NoSharing : public Sharing
{
  int getOwner(MeshEntity*) {return PCU_Comm_Self();} 
  bool isOwned(MeshEntity*) {return true;}
  virtual void getCopies(MeshEntity*, CopyArray&) {}
  bool isShared(MeshEntity*) {return false;}
};

Numbering* numberNodes(
    Mesh* mesh,
    const char* name,
    FieldShape* s,
    Sharing* shr, 
    bool delete_shr)
{
  Numbering* n = createNumbering(mesh,name,s,1);
  int i=0;
  for (int d=0; d < 4; ++d)
  {
    if ( ! s->hasNodesIn(d))
      continue;
    MeshIterator* it = mesh->begin(d);
    MeshEntity* e;
    while ((e = mesh->iterate(it)))
    {
      if ( ! shr->isOwned(e))
        continue;
      int nnodes = n->countNodesOn(e);
      for (int node=0; node < nnodes; ++node)
        number(n,e,node,0,i++);
    }
    mesh->end(it);
  }
  if (delete_shr) delete shr;
  return n;
}

Numbering* numberOwnedDimension(Mesh* mesh, const char* name, int dim,
    Sharing* shr)
{
  bool delete_shr=false;
  if (!shr) 
  {
    shr = getSharing(mesh);
    delete_shr=true;
  }
  return numberNodes(mesh, name, getConstant(dim), shr, delete_shr);
}

Numbering* numberOverlapDimension(Mesh* mesh, const char* name, int dim)
{
  FieldShape* s = getConstant(dim);
  Sharing* shr = new NoSharing();
  return numberNodes(mesh, name, s, shr, true);
}

Numbering* numberElements(Mesh* mesh, const char* name)
{
  return numberOwnedDimension(mesh, name, mesh->getDimension());
}

Numbering* numberOverlapNodes(Mesh* mesh, const char* name, FieldShape* s)
{
  if (!s)
    s = mesh->getShape();
  Sharing* shr = new NoSharing();
  return numberNodes(mesh, name, s, shr, true);
}

Numbering* numberOwnedNodes(
    Mesh* mesh,
    const char* name,
    FieldShape* s,
    Sharing* shr)
{
  if (!s)
    s = mesh->getShape();
  bool delete_shr=false;
  if (!shr)
  {
    shr = getSharing(mesh);
    delete_shr=true;
  }
  return numberNodes(mesh, name, s, shr, delete_shr);
}

class Counter : public FieldOp
{
  public:
    int count;
    FieldBase* field;
    virtual bool inEntity(MeshEntity* e)
    {
      if (field->getData()->hasEntity(e))
        count += field->countNodesOn(e);
      return false;
    }
    void run(FieldBase* f)
    {
      count = 0;
      field = f;
      apply(f);
    }
};

static int countFieldNodes(FieldBase* f)
{
  Counter counter;
  counter.run(f);
  return counter.count;
}

int countNodes(Numbering* n)
{
  return countFieldNodes(n);
}

int countNodes(GlobalNumbering* n)
{
  return countFieldNodes(n);
}

static void getFieldNodes(FieldBase* f, DynamicArray<Node>& nodes)
{
  Mesh* mesh = f->getMesh();
  FieldShape* s = f->getShape();
  nodes.setSize(countFieldNodes(f));
  size_t i = 0;
  for (int d=0; d < 4; ++d)
  {
    if ( ! s->hasNodesIn(d))
      continue;
    MeshIterator* it = mesh->begin(d);
    MeshEntity* e;
    while ((e = mesh->iterate(it)))
    {
      if ( ! f->getData()->hasEntity(e))
        continue;
      int nnodes = f->countNodesOn(e);
      for (int node=0; node < nnodes; ++node)
        nodes[i++] = Node(e,node);
    }
    mesh->end(it);
  }
  PCU_ALWAYS_ASSERT(i == nodes.getSize());
}

void getNodes(Numbering* n, DynamicArray<Node>& nodes)
{
  getFieldNodes(n,nodes);
}

typedef std::set<MeshEntity*> EntitySet;

static void getClosureEntitiesWithNodes(
    Mesh* m,
    MeshEntity* e,
    EntitySet& out,
    apf::FieldShape* s)
{
  int D = getDimension(m, e);
  for (int d=0; d <= D; ++d)
    if (s->hasNodesIn(d))
    {
      Downward de;
      int nde = m->getDownward(e,d,de);
      for (int i=0; i < nde; ++i)
        out.insert(de[i]);
    }
}

static void synchronizeEntitySet(
    Mesh* m,
    EntitySet& set)
{
  PCU_Comm_Begin();
  APF_ITERATE(EntitySet,set,it)
    if (m->isShared(*it))
    {
      Copies remotes;
      m->getRemotes(*it,remotes);
      APF_ITERATE(Copies,remotes,rit)
        PCU_COMM_PACK(rit->first,rit->second);
    }
  PCU_Comm_Send();
  while (PCU_Comm_Receive())
  {
    MeshEntity* e;
    PCU_COMM_UNPACK(e);
    set.insert(e);
  }
}

static void getNodesOnEntitySet(
    Mesh* m,
    EntitySet& s,
    DynamicArray<Node>& n,
    FieldShape* sh)
{
  size_t size = 0;
  APF_ITERATE(EntitySet,s,it) {
    int nen = sh->countNodesOn(m->getType(*it));
    size += nen;
  }
  n.setSize(size);
  size_t i=0;
  APF_ITERATE(EntitySet,s,it) {
    int nen = sh->countNodesOn(m->getType(*it));
    for (int j=0; j < nen; ++j)
      n[i++] = Node(*it,j);
  }
  PCU_ALWAYS_ASSERT(i==size);
}

void getNodesOnClosure(
    Mesh* m,
    ModelEntity* me,
    DynamicArray<Node>& on,
    FieldShape* sh)
{
  if (!sh)
    sh = m->getShape();
  int d = m->getModelType(me);
  MeshIterator* it = m->begin(d);
  MeshEntity* e;
  EntitySet s;
  while ((e = m->iterate(it)))
    if (m->toModel(e) == me)
      getClosureEntitiesWithNodes(m, e, s, sh);
  m->end(it);
  synchronizeEntitySet(m, s);
  getNodesOnEntitySet(m, s, on, sh);
}

GlobalNumbering* createGlobalNumbering(
    Mesh* mesh,
    const char* name,
    FieldShape* shape,
    int components)
{
  GlobalNumbering* n = new GlobalNumbering();
  n->init(name,mesh,shape,components);
  mesh->addGlobalNumbering(n);
  return n;
}

FieldShape* getShape(GlobalNumbering* n)
{
  return n->getShape();
}

const char* getName(GlobalNumbering* n)
{
  return n->getName();
}

Mesh* getMesh(GlobalNumbering* n)
{
  return n->getMesh();
}

int countComponents(GlobalNumbering* n)
{
  return n->countComponents();
}

void number(GlobalNumbering* n, Node node, long number, int component)
{
  n->set(node.entity,node.node,component,number);
}

void number(GlobalNumbering* n, MeshEntity* e, int node, long number)
{
  n->set(e, node, 0, number);
}

long getNumber(GlobalNumbering* n, Node node, int component)
{
  return n->get(node.entity,node.node,component);
}

long getNumber(GlobalNumbering* n, MeshEntity* e, int node, int component)
{
  return n->get(e,node,component);
}

int getElementNumbers(GlobalNumbering* n, MeshEntity* e,
    NewArray<long>& numbers)
{
  return n->getData()->getElementData(e,numbers);
}

static long exscan(long x)
{
  return PCU_Exscan_Long(x);
}

template <class T>
class Globalizer : public FieldOp
{
  public:
    T start;
    NumberingOf<T>* numbering;
    FieldDataOf<T>* data;
    std::vector<T> nodeData;
    virtual bool inEntity(MeshEntity* e)
    {
      if (data->hasEntity(e))
      {
        nodeData.resize(numbering->countNodesOn(e));
        data->get(e,&(nodeData[0]));
        for (size_t i=0; i < nodeData.size(); ++i)
          nodeData[i] += start;
        data->set(e,&(nodeData[0]));
      }
      return false;
    }
    void run(NumberingOf<T>* n)
    {
      numbering = n;
      data = n->getData();
      start = countFieldNodes(n);
      start = exscan(start);
      apply(n);
    }
};

static void globalize(GlobalNumbering* n)
{
  Globalizer<long> g;
  g.run(n);
}

void globalize(Numbering* n)
{
  Globalizer<int> g;
  g.run(n);
}

GlobalNumbering* makeGlobal(Numbering* n, bool destroy)
{
  std::string name = n->getName();
  name += "_global";
  Mesh* m = n->getMesh();
  FieldShape* s = n->getShape();
  GlobalNumbering* gn = createGlobalNumbering(
      m,name.c_str(),s);
  FieldDataOf<int>* nd = n->getData();
  for (int d = 0; d <= 3; ++d)
  {
    if (s->hasNodesIn(d))
    {
      MeshIterator* it = m->begin(d);
      MeshEntity* e;
      while ((e = m->iterate(it)))
        if (nd->hasEntity(e))
          for (int i = 0; i < n->countNodesOn(e); ++i)
          {
            apf::Node node(e,i);
            number(gn,node,getNumber(n,e,i,0));
          }
      m->end(it);
    }
  }
  if (destroy)
    apf::destroyNumbering(n);
  globalize(gn);
  return gn;
}

void synchronize(GlobalNumbering* n, Sharing* shr)
{
  synchronizeFieldData<long>(n->getData(), shr);
}

void destroyGlobalNumbering(GlobalNumbering* n)
{
  n->getMesh()->removeGlobalNumbering(n);
  delete n;
}

void getNodes(GlobalNumbering* n, DynamicArray<Node>& nodes)
{
  getFieldNodes(n,nodes);
}

}

