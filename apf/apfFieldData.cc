#include "apfFieldData.h"
#include "apfShape.h"
#include <PCU.h>

namespace apf {

FieldData::~FieldData()
{
}

template <class T>
void synchronizeFieldData(FieldDataOf<T>* data)
{
  FieldBase* f = data->getField();
  Mesh* m = f->getMesh();
  FieldShape* s = f->getShape();
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
          ( ! hasCopies(m,e))||
          ( ! isOriginal(m,e)))
        continue;
      int n = f->countValuesOn(e);
      NewArray<T> values(n);
      data->get(e,&(values[0]));
      Copies remotes;
      m->getRemotes(e,remotes);
      APF_ITERATE(Copies,remotes,it)
      {
        PCU_COMM_PACK(it->first,it->second);
        PCU_Comm_Pack(it->first,&(values[0]),n*sizeof(T));
      }
    }
    m->end(it);
    PCU_Comm_Send();
    while (PCU_Comm_Listen())
      while ( ! PCU_Comm_Unpacked())
      {
        MeshEntity* e;
        PCU_COMM_UNPACK(e);
        int n = f->countValuesOn(e);
        NewArray<T> values(n);
        PCU_Comm_Unpack(&(values[0]),n*sizeof(T));
        data->set(e,&(values[0]));
      }
  }
}

/* instantiate here */
template void synchronizeFieldData<int>(FieldDataOf<int>* data);
template void synchronizeFieldData<double>(FieldDataOf<double>* data);
template void synchronizeFieldData<long>(FieldDataOf<long>* data);

template <class T>
class Copy : public FieldOp
{
  public:
    Copy(FieldDataOf<T>* ld,
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
  Copy<T> copier(from,to);
  copier.run();
}

/* instantiate here */
template void copyFieldData<int>(
    FieldDataOf<int>* to, FieldDataOf<int>* from);
template void copyFieldData<double>(
    FieldDataOf<double>* to, FieldDataOf<double>* from);
template void copyFieldData<long>(
    FieldDataOf<long>* to, FieldDataOf<long>* from);

void accumulateFieldData(FieldDataOf<double>* data)
{
  FieldBase* f = data->getField();
  Mesh* m = f->getMesh();
  FieldShape* s = f->getShape();
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
          ( ! hasCopies(m,e))||
          (isOriginal(m,e)))
        continue;
      int n = f->countValuesOn(e);
      NewArray<double> values(n);
      data->get(e,&(values[0]));
      Copies remotes;
      m->getRemotes(e,remotes);
      APF_ITERATE(Copies,remotes,it)
      {
        PCU_COMM_PACK(it->first,it->second);
        PCU_Comm_Pack(it->first,&(values[0]),n*sizeof(double));
      }
    }
    m->end(it);
    PCU_Comm_Send();
    while (PCU_Comm_Listen())
      while ( ! PCU_Comm_Unpacked())
      {
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
  }
}

}
