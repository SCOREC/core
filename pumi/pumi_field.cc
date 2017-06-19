/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include "apf.h"
#include "apfShape.h"
#include "apfNumbering.h"
#include <pcu_util.h>
#include <PCU.h>
#include <cstdlib> // for malloc and free

//************************************
// Field shape and nodes
//************************************

pShape pumi_mesh_getShape (pMesh m)
{
  return m->getShape();
}

void pumi_mesh_setShape (pMesh m, pShape s, bool project)
{
  m->changeShape(s, project);
}

int pumi_shape_getNumNode (pShape s, int topo)
{
  return s->countNodesOn(topo);
}

bool pumi_shape_hasNode (pShape s, int topo)
{
  if (s->countNodesOn(topo)>0) 
    return true;
  return false;
}

void pumi_node_getCoord(pMeshEnt e, int i, double* xyz)
{
  Vector3 coord;
  pumi::instance()->mesh->getPoint(e, i, coord);
  for (int i=0; i<3; ++i)
    xyz[i] = coord[i]; 
}

void pumi_node_getCoordVector(pMeshEnt e, int i, Vector3& xyz)
{
  pumi::instance()->mesh->getPoint(e, i, xyz);
}


void pumi_node_setCoord(pMeshEnt e,int i, double* xyz)
{
  Vector3 coord;
  for (int k=0; k<3; ++k)
    coord[k] = xyz[k];
  pumi::instance()->mesh->setPoint(e, i, coord);
}

void pumi_node_setCoordVector(pMeshEnt e, int i, Vector3 const& xyz)
{
  pumi::instance()->mesh->setPoint(e, i, xyz);
}

/** \brief 3D vector cross product */
Vector3 pumi_vector3_cross(Vector3 const& a, Vector3 const& b)
{
  return apf::cross(a, b);
}


pShape pumi_shape_getLagrange (int order) { return apf::getLagrange(order); }
pShape pumi_shape_getSerendipity () { return apf::getSerendipity(); }
pShape pumi_shape_getConstant (int type) { return apf::getConstant(type); }
pShape pumi_shape_getIP (int dimension, int order) { return apf::getIPShape(dimension, order); }
pShape pumi_shape_getVoronoi (int dimension, int order) { return apf::getVoronoiShape(dimension, order); }
pShape pumi_shape_getIPFit(int dimension, int order) { return apf::getIPFitShape(dimension, order); }
pShape pumi_shape_getHierarchic (int order) { return apf::getHierarchic(order); }

//************************************
// Node numbering
//************************************

pGlobalNumbering pumi_numbering_createGlobal(pMesh m, const char* name, pShape shape, int num_component)
{
  if (!shape) shape= m->getShape();
  return apf::createGlobalNumbering(m, name, shape, num_component);
}

void pumi_numbering_deleteGlobal(pGlobalNumbering gn)
{
  apf::destroyGlobalNumbering(gn);
}

int pumi_mesh_getNumGlobalNumbering (pMesh m)
{
  return m->countGlobalNumberings();
}

pGlobalNumbering pumi_mesh_getGlobalNumbering (pMesh m, int i)
{
  return m->getGlobalNumbering(i);
}

void pumi_ment_setGlobalNumber(pMeshEnt e, pGlobalNumbering gn,
    int node, int component, long number)
{
  apf::Node n(e,node);
  apf::number(gn, n, component, number);
}

long pumi_ment_getGlobalNumber(pMeshEnt e, pGlobalNumbering gn, int node, int component)
{
  return apf::getNumber(gn, e, node, component);
}

pNumbering pumi_numbering_createLocalNode (pMesh m, const char* name, pShape shape)
{
  if (!shape) shape= m->getShape();
  return numberOverlapNodes(m, name, shape);
}

pNumbering pumi_numbering_createOwned (pMesh m, const char* name, int dim)
{
  return numberOwnedDimension(m, name, dim);
}

pNumbering pumi_numbering_create
   (pMesh m, const char* name, pShape shape, int num_component)
{
  if (!shape) shape= m->getShape();
  return createNumbering(m, name, shape, num_component);
}

pNumbering pumi_numbering_createOwnedNode (pMesh m, const char* name, pShape shape)
{
   if (!shape) shape= m->getShape();
   return numberOwnedNodes(m, name, shape);
}
void pumi_numbering_delete(pNumbering n)
{
  destroyNumbering(n);
}

int pumi_numbering_getNumNode(pNumbering n)
{
  return apf::countNodes(n);
}

//************************************
//  Field Management
//************************************ 

pField pumi_field_create(pMesh m, const char* name, int num_dof_per_ent, int type, pShape s)
{
  if (type==PUMI_PACKED)
    return apf::createPackedField(m, name, num_dof_per_ent, s);
  else
    return createGeneralField(m, name, type, num_dof_per_ent, s);
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
pShape pumi_field_getShape (pField f)
{
  return apf::getShape(f);
}

void pumi_field_delete(pField f)
{
  apf::destroyField(f);
}

void pumi_field_synchronize(pField f)
{  
  apf::synchronize(f, getSharing(getMesh(f)));
}

void pumi_field_accumulate(pField f)
{  
  apf::accumulate(f, getSharing(getMesh(f)));
}

void pumi_field_freeze(pField f)
{  
  if (!isFrozen(f))
    apf::freeze(f); 
}

void pumi_field_unfreeze(pField f)
{  
  if (isFrozen(f))
    apf::unfreeze(f); 
}


pField pumi_mesh_findField(pMesh m, const char* name)
{
  return m->findField(name);
}

int pumi_mesh_getNumField(pMesh m)
{
  return m->countFields();
}

pField pumi_mesh_getField(pMesh m, int i)
{
  return m->getField(i);
}

//*******************************************************
void pumi_ment_getField (pMeshEnt e, pField f, int i, double* dof_data)
//*******************************************************
{
  apf::getComponents(f, e, i, dof_data);
}    

void pumi_ment_setField (pMeshEnt e, pField f, int i, double* dof_data)
{
  apf::setComponents(f, e, i, dof_data);
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
    pMeshEnt e;
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

// VERIFY FIELDS
static void sendFieldData(pMesh m, pMeshEnt e, pField f, int nf)
{
  void* msg_send;
  pMeshEnt* s_ent;
  size_t msg_size;

  apf::Sharing* shr = getSharing(m);
    
  apf::FieldDataOf<double>* data = static_cast<apf::FieldDataOf<double>*>(f->getData());

  if ((!data->hasEntity(e))|| (!shr->isOwned(e)))
    return;

  int n = f->countValuesOn(e);

    msg_size=sizeof(pMeshEnt)+sizeof(int)+n*sizeof(double);
    apf::NewArray<double> values(n);
    data->get(e,&(values[0]));
    apf::CopyArray copies;
    shr->getCopies(e, copies);
    for (size_t i = 0; i < copies.getSize(); ++i)
    {
      int to = copies[i].peer;
      msg_send = malloc(msg_size);
      s_ent = (pMeshEnt*)msg_send; 
      *s_ent = copies[i].entity; 
      int *s_fieldid = (int*)((char*)msg_send + sizeof(pMeshEnt));
      s_fieldid[0] =  nf;
      double *s_data = (double*)((char*)msg_send+sizeof(pMeshEnt)+sizeof(int));
      for (int pos=0; pos<n; ++pos)
        s_data[pos]=values[pos];
      PCU_Comm_Write(to, (void*)msg_send, msg_size);
      free(msg_send);
    }

    if (m->isGhosted(e))
    {
      Copies g;
      m->getGhosts(e, g);
      APF_ITERATE(Copies, g, it)
      {
      int to = it->first;
      msg_send = malloc(msg_size);
      s_ent = (pMeshEnt*)msg_send; 
      *s_ent = it->second; 
      int *s_fieldid = (int*)((char*)msg_send + sizeof(pMeshEnt));
      s_fieldid[0] =  nf;
      double *s_data = (double*)((char*)msg_send+sizeof(pMeshEnt)+sizeof(int));
      for (int pos=0; pos<n; ++pos)
        s_data[pos]=values[pos];
      PCU_Comm_Write(to, (void*)msg_send, msg_size);
      free(msg_send);
      }
    } //if (m->isGhosted(e))
}

static void receiveFieldData(std::vector<pField>& fields, std::set<pField>& mismatch_fields)
{
  pField f;
  void *msg_recv;
  int pid_from;
  size_t msg_size;
  pMeshEnt e;

  while(PCU_Comm_Read(&pid_from, &msg_recv, &msg_size))
  {
    e = *((pMeshEnt*)msg_recv); 
    int *nf = (int*)((char*)msg_recv+sizeof(pMeshEnt)); 
    f = fields[*nf];

    int n = f->countValuesOn(e);

    double* r_values = (double*)((char*)msg_recv+sizeof(pMeshEnt)+sizeof(int)); 
    int num_data = (msg_size-sizeof(pMeshEnt)-sizeof(int))/sizeof(double);
    PCU_ALWAYS_ASSERT(n==num_data);
 
    if (mismatch_fields.find(f)!=mismatch_fields.end()) continue;

    apf::NewArray<double> values(n);
    f->getData()->get(e,&(values[0]));

    for (int i=0; i<n; ++i)
    {
      if (values[i]!=r_values[i])
      {
        mismatch_fields.insert(f);
        break;
      }
    }
  } // while
}

void pumi_field_verify(pMesh m, pField f)
{
  int n = m->countFields();
  if (!n) return;
  std::vector<pField> fields;
  if (f==NULL)
  {
    for (int i=0; i<n; ++i)
      fields.push_back(m->getField(i));
  }
  else 
    fields.push_back(f);

  if (!pumi_rank()) // master
  {
    printf("  - verifying fields: ");
    for (size_t nf = 0; nf < fields.size(); ++nf)
    {
      printf("%s", getName(fields[nf]));
      if (nf<fields.size()-1) printf(", ");      
    }
    printf("\n");
  }

  std::set<pField> mismatch_fields;
  for (size_t nf = 0; nf < fields.size(); ++nf)
  {
    for (int d=0; d < 4; ++d)
    {
      if (!static_cast<apf::FieldBase*>(fields[nf])->getShape()->hasNodesIn(d))
        continue;

      PCU_Comm_Begin();
      pMeshIter it = m->begin(d);
      pMeshEnt e;
      while ((e = m->iterate(it)))
        sendFieldData(m, e, f, nf);
      m->end(it);
      PCU_Comm_Send();
      receiveFieldData(fields,mismatch_fields); 
    }
  }
  int global_size = PCU_Max_Int((int)mismatch_fields.size());
  if (global_size&&!PCU_Comm_Self())
    for (std::set<pField>::iterator it=mismatch_fields.begin(); it!=mismatch_fields.end(); ++it)
      printf("%s: \"%s\" data mismatch over remote/ghost copies\n", __func__, getName(*it));
}

