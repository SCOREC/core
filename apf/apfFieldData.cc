#include <PCU.h>
#include "apfFieldData.h"
#include "apfShape.h"

namespace apf {

FieldData::~FieldData()
{
}

FieldData* FieldData::clone()
{
  abort();
}

void FieldData::rename(const char*)
{
  abort();
}

template <class T>
void synchronizeFieldData(FieldDataOf<T>* data, Sharing* shr)
{
  FieldBase* f = data->getField();
  Mesh* m = f->getMesh();
  FieldShape* s = f->getShape();
  if (!shr)
    shr = getSharing(m);
  for (int d=0; d < 4; ++d)
  {
    if ( ! s->hasNodesIn(d))
      continue;
    MeshEntity* e;
    MeshIterator* it = m->begin(d);
    PCU_Comm_Begin();
    while ((e = m->iterate(it)))
    {
      if (( ! data->hasEntity(e))||
          ( ! shr->isOwned(e)))
        continue;
      int n = f->countValuesOn(e);
      NewArray<T> values(n);
      data->get(e,&(values[0]));
      CopyArray copies;
      shr->getCopies(e, copies);
      for (size_t i = 0; i < copies.getSize(); ++i)
      {
        PCU_COMM_PACK(copies[i].peer, copies[i].entity);
        PCU_Comm_Pack(copies[i].peer, &(values[0]), n*sizeof(T));
      }
    }
    m->end(it);
    PCU_Comm_Send();
    while (PCU_Comm_Receive())
    {
      MeshEntity* e;
      PCU_COMM_UNPACK(e);
      int n = f->countValuesOn(e);
      NewArray<T> values(n);
      PCU_Comm_Unpack(&(values[0]),n*sizeof(T));
      data->set(e,&(values[0]));
    }
  }
  delete shr;
}

/* instantiate here */
template void synchronizeFieldData<int>(FieldDataOf<int>*, Sharing*);
template void synchronizeFieldData<double>(FieldDataOf<double>*, Sharing*);
template void synchronizeFieldData<long>(FieldDataOf<long>*, Sharing*);

template <class T>
class CopyOp : public FieldOp
{
  public:
    CopyOp(FieldDataOf<T>* ld,
         FieldDataOf<T>* rd)
    {
      from = ld;
      to = rd;
    }
    bool inEntity(MeshEntity* e)
    {
      if (from->hasEntity(e))
      {
        int n = from->getField()->countValuesOn(e);
        NewArray<T> v(n);
        from->get(e,&(v[0]));
        to->set(e,&(v[0]));
      }
      return false;
    }
    void run() {apply(to->getField());}
    FieldDataOf<T>* from;
    FieldDataOf<T>* to;
};

template <class T>
void copyFieldData(FieldDataOf<T>* to, FieldDataOf<T>* from)
{
  CopyOp<T> copier(from,to);
  copier.run();
}

/* instantiate here */
template void copyFieldData<int>(
    FieldDataOf<int>* to, FieldDataOf<int>* from);
template void copyFieldData<double>(
    FieldDataOf<double>* to, FieldDataOf<double>* from);
template void copyFieldData<long>(
    FieldDataOf<long>* to, FieldDataOf<long>* from);

void accumulateFieldData(FieldDataOf<double>* data, Sharing* shr)
{
  FieldBase* f = data->getField();
  Mesh* m = f->getMesh();
  FieldShape* s = f->getShape();
  if (!shr)
    shr = getSharing(m);
  for (int d=0; d < 4; ++d)
  {
    if ( ! s->hasNodesIn(d))
      continue;
    MeshEntity* e;
    MeshIterator* it = m->begin(d);
    PCU_Comm_Begin();
    while ((e = m->iterate(it)))
    {
      if (( ! data->hasEntity(e)) ||
          (shr->isOwned(e)))
        continue; /* non-owners send to owners */
      CopyArray copies;
      shr->getCopies(e, copies);
      int n = f->countValuesOn(e);
      NewArray<double> values(n);
      data->get(e,&(values[0]));
      /* actually, non-owners send to all others,
         since apf::Sharing doesn't identify the owner */
      for (size_t i = 0; i < copies.getSize(); ++i)
      {
        PCU_COMM_PACK(copies[i].peer, copies[i].entity);
        PCU_Comm_Pack(copies[i].peer, &(values[0]), n*sizeof(double));
      }
    }
    m->end(it);
    PCU_Comm_Send();
    while (PCU_Comm_Listen())
      while ( ! PCU_Comm_Unpacked())
      { /* receive and add. we only care about correctness
           on the owners */
        MeshEntity* e;
        PCU_COMM_UNPACK(e);
        int n = f->countValuesOn(e);
        NewArray<double> values(n);
        NewArray<double> inValues(n);
        PCU_Comm_Unpack(&(inValues[0]),n*sizeof(double));
        data->get(e,&(values[0]));
        for (int i = 0; i < n; ++i)
          values[i] += inValues[i];
        data->set(e,&(values[0]));
      }
  } /* broadcast back out to non-owners */
  synchronizeFieldData(data, shr);
}

template <class T>
void FieldDataOf<T>::setNodeComponents(MeshEntity* e, int node,
    T const* components)
{
  int n = field->countNodesOn(e);
  if (n==1)
    return set(e,components);
  int nc = field->countComponents();
  NewArray<T> allComponents(nc*n);
  if (this->hasEntity(e))
    get(e,&(allComponents[0]));
  for (int i=0; i < nc; ++i)
    allComponents[node*nc+i] = components[i];
  set(e,&(allComponents[0]));
}

template <class T>
void FieldDataOf<T>::getNodeComponents(MeshEntity* e, int node, T* components)
{
  int n = field->countNodesOn(e);
  if (n==1)
    return get(e,components);
  int nc = field->countComponents();
  NewArray<T> allComponents(nc*n);
  get(e,&(allComponents[0]));
  for (int i=0; i < nc; ++i)
    components[i] = allComponents[node*nc+i];
}

template <class T>
void reorderData(T const dataIn[], T dataOut[], int const order[], int nc, int nn)
{
  for (int i = 0; i < nn; ++i) {
    int oi = order[i];
    for (int j = 0; j < nc; ++j)
      dataOut[oi * nc + j] = dataIn[i * nc + j];
  }
}

template <class T>
int FieldDataOf<T>::getElementData(MeshEntity* entity, NewArray<T>& data)
{
  Mesh* mesh = field->getMesh();
  int t = mesh->getType(entity);
  int ed = Mesh::typeDimension[t];
  FieldShape* fs = field->getShape();
  EntityShape* es = fs->getEntityShape(t);
  int nc = field->countComponents();
  int nen = es->countNodes();
  data.allocate(nc * nen);
  /* try to minimize the amount of
     reallocation of these two.
     when there is no reordering, they
     shouln't allocate at all */
  apf::DynamicArray<int> order;
  apf::DynamicArray<T> adata;
  int n = 0;
  for (int d = 0; d <= ed; ++d)
  {
    if (fs->hasNodesIn(d))
    {
      Downward a;
      int na = mesh->getDownward(entity,d,a);
      for (int i = 0; i < na; ++i)
      {
        int nan = fs->countNodesOn(mesh->getType(a[i]));
        if (nan > 1 && ed != d) {
          order.setSize(nen); /* nen >= nan */
          adata.setSize(nen); /* setSize is no-op for the same size */
          es->alignSharedNodes(mesh, entity, a[i], &order[0]);
          get(a[i], &adata[0]);
          reorderData<T>(&adata[0], &data[n], &order[0], nc, nan);
        } else {
          get(a[i], &data[n]);
        }
        n += nc * nan;
      }
    }
  }
  assert(n == nc * nen);
  return n;
}

template class FieldDataOf<double>;
template class FieldDataOf<int>;
template class FieldDataOf<long>;

}
