#include <PCU.h>
#include "apfFieldData.h"
#include "apfShape.h"
#include <pcu_util.h>
#include <cstdlib>
#include <iostream>

namespace apf {

FieldData::~FieldData()
{
}

void FieldData::rename(const char*)
{
  abort();
}

template <class T>
void synchronizeFieldData(FieldDataOf<T>* data, Sharing* shr, bool delete_shr)
{
  FieldBase* f = data->getField();
  Mesh* m = f->getMesh();
  FieldShape* s = f->getShape();
  if (!shr)
  {
    shr = getSharing(m);
    delete_shr=true;
  }
  for (int d=0; d < 4; ++d)
  {
    if ( ! s->hasNodesIn(d))
      continue;
    MeshEntity* e;
    MeshIterator* it = m->begin(d);
    m->getPCU()->Begin();
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
        m->getPCU()->Pack(copies[i].peer, copies[i].entity);
        m->getPCU()->Pack(copies[i].peer, &(values[0]), n*sizeof(T));
      }
      apf::Copies ghosts;  
      if (m->getGhosts(e, ghosts))
      APF_ITERATE(Copies, ghosts, it)
      {
        m->getPCU()->Pack(it->first, it->second);
        m->getPCU()->Pack(it->first, &(values[0]), n*sizeof(T));
      }
    }
    m->end(it);
    m->getPCU()->Send();
    while (m->getPCU()->Receive())
    {
      MeshEntity* e;
      PCU_COMM_UNPACK(e);
      int n = f->countValuesOn(e);
      NewArray<T> values(n);
      m->getPCU()->Unpack(&(values[0]),n*sizeof(T));
      data->set(e,&(values[0]));
    }
  }
  if (delete_shr) delete shr;
}

/* instantiate here */
template void synchronizeFieldData<int>(FieldDataOf<int>*, Sharing*, bool);
template void synchronizeFieldData<double>(FieldDataOf<double>*, Sharing*, bool);
template void synchronizeFieldData<long>(FieldDataOf<long>*, Sharing*, bool);

void reduceFieldData(FieldDataOf<double>* data, Sharing* shr, bool delete_shr, const ReductionOp<double>& reduce_op /* =ReductionSum<double>() */)
{
  FieldBase* f = data->getField();
  Mesh* m = f->getMesh();
  FieldShape* s = f->getShape();
  if (!shr)
  {
    shr = getSharing(m);
    delete_shr=true;
  }
  for (int d=0; d < 4; ++d)
  {
    if ( ! s->hasNodesIn(d))
      continue;

    MeshEntity* e;
    MeshIterator* it = m->begin(d);
    m->getPCU()->Begin();
    while ((e = m->iterate(it)))
    {
      /* send to all parts that can see this entity */
      if ( ! data->hasEntity(e) )
        continue;
 
      if (m->isGhost(e) && shr->isShared(e))
      {
        // zero out ghost values (because we reduce only over non-ghost values)
        int n = f->countValuesOn(e);
        NewArray<double> values(n);
        for (int i=0; i < n; ++i)
          values[i] = reduce_op.getNeutralElement();

        data->set(e, &(values[0]));
        continue;
      }

      // copies
      CopyArray copies;
      shr->getCopies(e, copies);
      int n = f->countValuesOn(e);
      NewArray<double> values(n);
      data->get(e,&(values[0]));

      for (size_t i = 0; i < copies.getSize(); ++i)
      {
        m->getPCU()->Pack(copies[i].peer, copies[i].entity);
        m->getPCU()->Pack(copies[i].peer, &(values[0]), n*sizeof(double));
      }

      // ghosts - only do them if this entity is on a partition boundary
      if (copies.getSize() > 0)
      {
        apf::Copies ghosts;
        if (m->getGhosts(e, ghosts))
        APF_ITERATE(Copies, ghosts, it2)
        {
          m->getPCU()->Pack(it2->first, it2->second);
          m->getPCU()->Pack(it2->first, &(values[0]), n*sizeof(double));
        }
      }
    }
    m->end(it);

    m->getPCU()->Send();
    while (m->getPCU()->Listen())
      while ( ! m->getPCU()->Unpacked())
      { /* receive and add. we only care about correctness
           on the owners */
        MeshEntity* e;
        PCU_COMM_UNPACK(e);
        int n = f->countValuesOn(e);
        NewArray<double> values(n);
        NewArray<double> inValues(n);
        m->getPCU()->Unpack(&(inValues[0]),n*sizeof(double));
        data->get(e,&(values[0]));
        for (int i = 0; i < n; ++i)
          values[i] = reduce_op.apply(values[i], inValues[i]);
        data->set(e,&(values[0]));
      }
  }

  // every partition did the reduction,s o no need to broadcast the result
  if (delete_shr) delete shr;
}

template <class T>
void FieldDataOf<T>::setNodeComponents(MeshEntity* e, int node,
    T const* components)
{
  int n = field->countNodesOn(e);
  if (n==1) {
    PCU_ALWAYS_ASSERT(node == 0);
    return set(e,components);
  }
  PCU_ALWAYS_ASSERT(node >= 0);
  PCU_ALWAYS_ASSERT(node < n);
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
  if (n==1) {
    PCU_ALWAYS_ASSERT(node == 0);
    return get(e,components);
  }
  PCU_ALWAYS_ASSERT(node >= 0);
  PCU_ALWAYS_ASSERT(node < n);
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

// This is only used to reorder the data for interior face nodes on a face of
// a Nedelec tet, where each node on the face contains 2 dof values.
template <class T>
void reorderDataNedelec(
    T const dataIn[],
    T dataOut[],
    int const order[],
    int nc,
    int nn,
    int type)
{
  if (type == Mesh::TRIANGLE)
    for (int i = 0; i < nn; ++i) {
      if(order[2*nn+i])
      {
	dataOut[i*nc] = (order[i] >= 0) ?
	  dataIn[ order[i] ] : -dataIn[ -(order[i]+1) ];
      }
      else
	dataOut[i*nc] = 0.;

      if(order[3*nn+i])
      {
	dataOut[i*nc] += (order[1*nn+i] >= 0) ?
	  dataIn[ order[1*nn+i] ] : -dataIn[ -(order[1*nn+i]+1) ];
      }
    }
  else if (type == Mesh::EDGE)
    for (int i = 0; i < nn; ++i) {
      int oi = order[i] >= 0 ? order[i] : -(order[i]+1);
      for (int j = 0; j < nc; ++j)
	dataOut[oi * nc + j] =
	  order[i] >= 0 ? dataIn[i * nc + j] : -dataIn[i * nc + j];
    }
  else
    PCU_ALWAYS_ASSERT_VERBOSE(0,
    	"type has to be Mesh::EDGE or Mesh::TRIANGLE!");
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
  for (int d = 0; d <= ed; ++d) {
    if (!fs->hasNodesIn(d)) continue;
    Downward a;
    int na = mesh->getDownward(entity,d,a);
    for (int i = 0; i < na; ++i) {
      int nan = fs->countNodesOn(mesh->getType(a[i]));
      // for vector shapes (i.e., nedelec) direction matters for nan>=1
      if (fs->isVectorShape()) {
        if (nan >= 1 && ed != d) {
          // The 1st nen ints is the 1st contribution
          // The 2nd nen ints is the 2nd contribution
          // The 3rd nen ints tells whether to add 1st contribution
          // The 4th nen ints tells whether to add 2st contribution
          order.setSize(4*nen);
          adata.setSize(nen);
          es->alignSharedNodes(mesh, entity, a[i], &order[0]);
          get(a[i], &adata[0]);
          // We would want to have a different reorder here to handle
          // the fact that order now includes some extra info
          int dtype = mesh->getType(a[i]);
          reorderDataNedelec<T>(&adata[0], &data[n], &order[0], nc, nan, dtype);
        }
        // this else is required to add the dofs associated with the tet
        else
          get(a[i], &data[n]);
      }
      // for non vector shapes direction matters for nan>1
      else {
	if (nan > 1 && ed != d) { /* multiple shared nodes, check alignment */
	  order.setSize(nen); /* nen >= nan */
	  adata.setSize(nen); /* setSize is no-op for the same size */
	  // Note: The above efficiency consideration does not account for the
	  // fact that nc might be very large (e.g. nc = 9 for matrix fields)
	  // and for such cases setting the size of adata to "nen" is not enough.
	  // Hence the need for the following line.
	  if (nan*nc > nen)
	    adata.setSize(nan*nc);
	  es->alignSharedNodes(mesh, entity, a[i], &order[0]);
	  get(a[i], &adata[0]);
	  reorderData<T>(&adata[0], &data[n], &order[0], nc, nan);
	} else if (nan) { /* non-zero set of nodes, either one
			    or not shared */
	  get(a[i], &data[n]);
	}
      }
      n += nc * nan;
    }
  }
  PCU_ALWAYS_ASSERT(n == nc * nen);
  return n;
}

template class FieldDataOf<double>;
template class FieldDataOf<int>;
template class FieldDataOf<long>;

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
void copyFieldData(FieldDataOf<T>* from, FieldDataOf<T>* to)
{
  CopyOp<T> copier(from,to);
  copier.run();
}

/* instantiate here */
template void copyFieldData<int>(
    FieldDataOf<int>* from, FieldDataOf<int>* to);
template void copyFieldData<double>(
    FieldDataOf<double>* from, FieldDataOf<double>* to);
template void copyFieldData<long>(
    FieldDataOf<long>* from, FieldDataOf<long>* to);

template <class T>
class MultiplyOp : public FieldOp
{
  public:
    MultiplyOp(FieldDataOf<T>* ld, T d,
         FieldDataOf<T>* rd)
    {
      from = ld;
      to = rd;
      mult = d;
    }
    bool inEntity(MeshEntity* e)
    {
      if (from->hasEntity(e))
      {
        int n = from->getField()->countValuesOn(e);
        NewArray<T> v(n);
        from->get(e,&(v[0]));
        for (int i=0; i<n; ++i)
          v[i]*=mult;
        to->set(e,&(v[0]));
      }
      return false;
    }
    void run() {apply(to->getField());}
    FieldDataOf<T>* from;
    T mult;
    FieldDataOf<T>* to;
};

template <class T>
void multiplyFieldData(FieldDataOf<T>* from, T d, FieldDataOf<T>* to)
{
  MultiplyOp<T> multiplier(from,d,to);
  multiplier.run();
}

/* instantiate here */
template void multiplyFieldData<int>(
    FieldDataOf<int>* from, int d, FieldDataOf<int>* to);
template void multiplyFieldData<double>(
    FieldDataOf<double>* from, double d, FieldDataOf<double>* to);
template void multiplyFieldData<long>(
    FieldDataOf<long>* from, long d, FieldDataOf<long>* to);


template <class T>
class AddOp : public FieldOp
{
  public:
    AddOp(FieldDataOf<T>* ld1, FieldDataOf<T>* ld2,
         FieldDataOf<T>* rd)
    {
      from1 = ld1;
      from2 = ld2;
      to = rd;
    }
    bool inEntity(MeshEntity* e)
    {
      if (from1->hasEntity(e) && from2->hasEntity(e))
      {
        int n = from1->getField()->countValuesOn(e);
        NewArray<T> v1(n);
        from1->get(e,&(v1[0]));
        NewArray<T> v2(n);
        from2->get(e,&(v2[0]));
        for (int i=0; i<n; ++i)
          v1[i]+=v2[i];
        to->set(e,&(v1[0]));
      }
      return false;
    }
    void run() {apply(to->getField());}
    FieldDataOf<T>* from1;
    FieldDataOf<T>* from2;
    FieldDataOf<T>* to;
};

template <class T>
void addFieldData(FieldDataOf<T>* from1, FieldDataOf<T>* from2, FieldDataOf<T>* to)
{
  AddOp<T> adder(from1,from2,to);
  adder.run();
}

/* instantiate here */
template void addFieldData<int>(
    FieldDataOf<int>*, FieldDataOf<int>*,FieldDataOf<int>*);
template void addFieldData<double>(
    FieldDataOf<double>*, FieldDataOf<double>*, FieldDataOf<double>*);
template void addFieldData<long>(
    FieldDataOf<long>*, FieldDataOf<long>*, FieldDataOf<long>*);

}
