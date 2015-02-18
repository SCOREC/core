#ifndef APFSIMDATAOF_H
#define APFSIMDATAOF_H

#include "apfFieldData.h"
#include <SimField.h>


namespace apf {

template <class T>
class SIMDataOf : public FieldDataOf<T>
{
public:
  virtual void init(FieldBase * f)
  {
    FieldData::field = f;
    mesh = f->getMesh();
    pf = PolyField_new(f->getShape()->getOrder(), 0);
    fd = Field_new(static_cast<MeshSIM*>(mesh)->mesh,
		   f->countComponents(),
		   "apf_field_data",
		   f->getName(),
		   ShpLagrange,
		   1, // surfcont (always 1 in current simModSuite)
		   1, // num_time_derivatives to keep available
		   1, // num_section
		   pf);
    Field_apply(fd,mesh->getDimension(),NULL);
  }
  virtual bool hasEntity(MeshEntity * e)
  {
    return true;
  }
  virtual void removeEntity(MeshEntity * e)
  {}
  virtual void get(MeshEntity * e, T * data)
  {
    pEntity ent = reinterpret_cast<pEntity>(e);
    pDofGroup dof;
    int num_read = 0;
    for(int eindex = 0; (dof = Field_entDof(fd,ent,eindex)) ; eindex++)
    {
      int dofs_per_node = DofGroup_numComp(dof);
      DofGroup_values(dof,0,&data[num_read]);
      num_read += dofs_per_node;
    }
  }
  virtual void set(MeshEntity * e, T const * data)
  {
    pEntity ent = reinterpret_cast<pEntity>(e);
    pDofGroup dof;
    int num_wrote = 0;
    for(int eindex = 0; (dof = Field_entDof(fd,ent,eindex)); eindex++)
    {
      int dofs_per_node = DofGroup_numComp(dof);
      for(int ii = 0; ii < dofs_per_node; ii++)
	DofGroup_setValue(dof,ii,0,data[num_wrote+ii]);
      num_wrote += dofs_per_node;
    }
  }
  virtual bool isFrozen() {return false;}
  virtual FieldData * clone()
  {
    return NULL;
  }

private:
  Mesh * mesh;
  pField fd;
  pPolyField pf;
};

} // namespace apf

#endif
