/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include "apf.h"
#include "apfShape.h"

void pumi_mesh_getField(pMesh m, std::vector<pField>& fields)
{
  for (int i=0; i<m->countFields(); ++i)
    fields.push_back(m->getField(i));
}

pField pumi_field_create(pMesh m, const char* name, int num_dof_per_ent, int type, pShape s)
{
  if (type==apf::PACKED)
    return apf::createPackedField(m, name, num_dof_per_ent, s);
  else
  {
    return createGeneralField(m, name, type, num_dof_per_ent, s);
  }
}

int pumi_field_getSize(pField f)
{  
  return apf::countComponents(f); 
}

int pumi_field_getType(pField f)
{ 
  return apf::getValueType(f); 
}

std::string pumi_field_getName(pField f)
{ 
  return apf::getName(f); 
}

void pumi_field_delete(pField f)
{
  apf::destroyField(f);
}

void pumi_field_synchronize(apf::Field* f)
{  
  apf::synchronize(f, getSharing(getMesh(f)));
}

void pumi_field_accumulate(apf::Field* f)
{  
  apf::accumulate(f, getSharing(getMesh(f)));
}

void pumi_field_freeze(apf::Field* f)
{  
  if (!isFrozen(f))
    apf::freeze(f); 
}

void pumi_field_unfreeze(apf::Field* f)
{  
  if (isFrozen(f))
    apf::unfreeze(f); 
}

//*******************************************************
void pumi_ment_getField (pMeshEnt e, pField f, double* dof_data)
//*******************************************************
{
  apf::getComponents(f, e, 0, dof_data);
}    

void pumi_ment_setField (pMeshEnt e, pField f, double* dof_data)
{
  apf::setComponents(f, e, 0, dof_data);
}

pField pumi_mesh_findField(pMesh m, const char* name)
{
  return m->findField(name);
}

#include "apfField.h"
#include "apfFieldData.h"
#include <iostream>
//*******************************************************
void pumi_field_print(pField f)
//*******************************************************
{
  pumi_sync();
  apf::Mesh*  m = getMesh(f);
  if (!m->findTag("global_id")) pumi_mesh_createGlobalID((pMesh)m);

  for (int d=0; d < 4; ++d)
  {
    if (!static_cast<apf::FieldBase*>(f)->getShape()->hasNodesIn(d))
      continue;

    apf::FieldDataOf<double>* data = static_cast<apf::FieldDataOf<double>*>(f->getData());
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it)))
    {
      if (!data->hasEntity(e))
        continue;

      int n = f->countValuesOn(e);
      apf::NewArray<double> dof_data(n);
      data->get(e,&(dof_data[0]));
      switch (n)
      {
        case 1: {
          std::cout<<"[p"<<pumi_rank()<<"] field "<<getName(f)
		     <<"/ent "<<e<<" id "<<pumi_ment_getGlobalID(e)
		     <<": ["<<dof_data[0]
		     <<"]\n";
        break;}
      case 2: {     
	std::cout<<"[p"<<pumi_rank()<<"] field "<<getName(f)
		     <<"/ent "<<e<<" id "<<pumi_ment_getGlobalID(e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<"]\n";
        break;}
      case 3: {
	std::cout<<"[p"<<pumi_rank()<<"] field "<<getName(f)
		     <<"/ent "<<e<<" id "<<pumi_ment_getGlobalID(e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<"]\n";
        break;}
    case 4: {
	std::cout<<"[p"<<pumi_rank()<<"] field "<<getName(f)
		     <<"/ent "<<e<<" id "<<pumi_ment_getGlobalID(e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<"]\n";
 
        break; }
      case 6: {
	std::cout<<"[p"<<pumi_rank()<<"] field "<<getName(f)
		     <<"/ent "<<e<<" id "<<pumi_ment_getGlobalID(e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<"]\n";
        break; }
      case 8: {
	std::cout<<"[p"<<pumi_rank()<<"] field "<<getName(f)
		     <<"/ent "<<e<<" id "<<pumi_ment_getGlobalID(e)
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<"]\n";
        break; }
      case 12: {
	std::cout<<"[p"<<pumi_rank()<<"] field "<<getName(f)
		     <<"/ent "<<e<<" id "<<pumi_ment_getGlobalID(e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<", "<<dof_data[8]
		     <<", "<<dof_data[9]
		     <<", "<<dof_data[10]
		     <<", "<<dof_data[11]
		     <<"]\n";
        break; }
      case 18: {
	std::cout<<"[p"<<pumi_rank()<<"] field "<<getName(f)
		     <<"/ent "<<e<<" id "<<pumi_ment_getGlobalID(e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<", "<<dof_data[8]
		     <<", "<<dof_data[9]
		     <<", "<<dof_data[10]
		     <<", "<<dof_data[11]
		     <<", "<<dof_data[12]
		     <<", "<<dof_data[13]
		     <<", "<<dof_data[14]
		     <<", "<<dof_data[15]
		     <<", "<<dof_data[16]
		     <<", "<<dof_data[17]
		     <<"]\n";
        break; }
      case 24: {
	std::cout<<"[p"<<pumi_rank()<<"] field "<<getName(f)
		     <<"/ent "<<e<<" id "<<pumi_ment_getGlobalID(e)
		     <<": ["<<dof_data[0]
		     <<", "<<dof_data[1]
		     <<", "<<dof_data[2]
		     <<", "<<dof_data[3]
		     <<", "<<dof_data[4]
		     <<", "<<dof_data[5]
		     <<", "<<dof_data[6]
		     <<", "<<dof_data[7]
		     <<", "<<dof_data[8]
		     <<", "<<dof_data[9]
		     <<", "<<dof_data[10]
		     <<", "<<dof_data[11]
		     <<", "<<dof_data[12]
		     <<", "<<dof_data[13]
		     <<", "<<dof_data[14]
		     <<", "<<dof_data[15]
		     <<", "<<dof_data[16]
		     <<", "<<dof_data[17]
		     <<", "<<dof_data[18]
		     <<", "<<dof_data[19]
		     <<", "<<dof_data[20]
		     <<", "<<dof_data[21]
		     <<", "<<dof_data[22]
		     <<", "<<dof_data[23]
		     <<"]\n";
        break; }
      default: if (!pumi_rank()) std::cout<<__func__<<" failed for field "
               <<getName(f)<<": does support "<<n<<" dofs\n";
               break;
      } // switch
    } // while
    m->end(it);
  }
}

