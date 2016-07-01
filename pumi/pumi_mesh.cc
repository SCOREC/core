/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#include "pumi.h"
#include <mpi.h>
#include <parma.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <apfZoltan.h>
#include <assert.h>
#include <iostream>
#include <string.h>

//*******************************************************
void destroy_node_global_numbering(apf::Mesh2* m)
//*******************************************************
{
  if (m->findField("node own partid field"))
    destroyField(m->findField("node own partid field"));
  if (m->findField("node global id field"))
    destroyField(m->findField("node global id field"));
}

//*******************************************************
void generate_node_global_numbering(apf::Mesh2* m)
//*******************************************************
{
//  if (!PCU_Comm_Self())  std::cout<<"[M3D-C1 INFO] ***** GENERATING GLOBAL NODE NUMBERING ***** \n"; 
  destroy_node_global_numbering(m);

  double id[1];
  double own_partid[1];
  apf::Field* node_ownpid_f = createPackedField(m, "node own partid field", 1);
  apf::freeze(node_ownpid_f);
  apf::Field* node_globalid_f = createPackedField(m, "node global id field", 1);
  apf::freeze(node_globalid_f);

  // count #own_vtx
  int num_own_ent=0;
  pMeshEnt e;
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    if (m->getOwner(e)==PCU_Comm_Self())
      ++num_own_ent;
  }
  m->end(it);

  // generate global node_id
  pumi::instance()->num_own_vtx=num_own_ent;
  PCU_Exscan_Ints(&num_own_ent,1);
  int start=num_own_ent;

  PCU_Comm_Begin();

  it = m->begin(0);
  while ((e = m->iterate(it)))
  {
    own_partid[0]=(double)m->getOwner(e); 
    setComponents(node_ownpid_f, e, 0, own_partid);    
    if ((int)(own_partid[0])!=PCU_Comm_Self()) continue;
    id[0] = (double) start;
    setComponents(node_globalid_f, e, 0, id);
    pCopies remotes;
    m->getRemotes(e,remotes);
    APF_ITERATE(pCopies,remotes,it)
    {
      PCU_COMM_PACK(it->first,it->second);
      PCU_Comm_Pack(it->first,&start,sizeof(int));
    }
    ++start;
  }
  m->end(it);
  PCU_Comm_Send();

  int value;
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      pMeshEnt r;
      PCU_COMM_UNPACK(r);
      PCU_Comm_Unpack(&value,sizeof(int));
      id[0] = (double) value;
      setComponents(node_globalid_f, r, 0, id);
    }
}


// *********************************************************
pumi::pumi(): mesh(NULL), model(NULL) {}

pumi::~pumi()
{
  delete _instance;
  _instance = NULL;
}

pumi* pumi::_instance=NULL;
pumi* pumi::instance()
{
  if (_instance==NULL)
    _instance = new pumi();
  return _instance;
}

apf::Migration* getPlan(apf::Mesh* m, int num_target_part)
{
  apf::Splitter* splitter = apf::makeZoltanSplitter(
      m, apf::GRAPH, apf::PARTITION, false);
  pTag weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.05, num_target_part);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;
  return plan;
}

void split_comm(int num_out_comm)
{
  int self = PCU_Comm_Self();
  int group_id = self % num_out_comm;
  int in_group_rank = self / num_out_comm;
  MPI_Comm groupComm;
  MPI_Comm_split(PCU_Get_Comm(), group_id, in_group_rank, &groupComm);
  PCU_Switch_Comm(groupComm);
}

void merge_comm(MPI_Comm oldComm)
{
  MPI_Comm prevComm = PCU_Get_Comm();
  PCU_Switch_Comm(oldComm);
  MPI_Comm_free(&prevComm);
}

// *********************************************************
void assign_uniq_partbdry_id(pMesh m, int dim, pTag partbdry_id_tag)
// *********************************************************
{
  pMeshEnt e;

  int own_partid, num_own_partbdry=0, myrank=PCU_Comm_Self();

  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
  {
    own_partid = m->getOwner(e);
    if (own_partid==myrank && m->isShared(e))
      ++num_own_partbdry;
  }
  m->end(it);

  PCU_Exscan_Ints(&num_own_partbdry,1);
  int initial_id=num_own_partbdry;

  PCU_Comm_Begin();
  it = m->begin(dim);
  while ((e = m->iterate(it)))
  {
    own_partid = m->getOwner(e);
    if (own_partid==myrank && m->isShared(e))
    {
      m->setIntTag(e, partbdry_id_tag, &initial_id);
      pCopies remotes;
      m->getRemotes(e, remotes);
      APF_ITERATE(pCopies, remotes, it)
      {
        PCU_COMM_PACK(it->first, it->second);
        PCU_Comm_Pack(it->first, &initial_id, sizeof(int));
      }
      ++initial_id;
    }
  }
  m->end(it);

  PCU_Comm_Send();
  int global_id;
  while (PCU_Comm_Listen())
    while (!PCU_Comm_Unpacked())
    {
      pMeshEnt remote_ent;
      PCU_COMM_UNPACK(remote_ent);
      PCU_Comm_Unpack(&global_id, sizeof(int));
      m->setIntTag(remote_ent, partbdry_id_tag, &global_id);
    }
}

void set_remote(pMesh m, pMeshEnt e, int p, pMeshEnt r)
{
  pCopies remotes;
  m->getRemotes(e,remotes);
  bool found=false;
  APF_ITERATE(pCopies, remotes, it)
    if (it->first==p) 
    {
      found=true;
      break;
    }
  if (found)
  {
    remotes[p] = r;
    m->setRemotes(e,remotes);
  }
  else
    m->addRemote(e, p, r);
}


// *********************************************************
void stitch_link(pMesh m, pTag partbdry_id_tag, 
                      std::map<int, pMeshEnt>* partbdry_entities)
// *********************************************************
{
  PCU_Comm_Begin();
  pMeshEnt e;
  
  int ent_info[3];
  int global_id; 
  for (int dim=0; dim<2; ++dim)
  {
    for (std::map<int, pMeshEnt>::iterator ent_it = partbdry_entities[dim].begin();
         ent_it!= partbdry_entities[dim].end(); ++ent_it)
    {    
      e = ent_it->second;
      pCopies remotes;
      m->getRemotes(e, remotes);
      m->getIntTag(e, partbdry_id_tag, &global_id);
      ent_info[0] = dim;
      ent_info[1] = global_id; 
      ent_info[2] = PCU_Comm_Self(); 
      APF_ITERATE(apf::Copies, remotes, it)
      {  
        PCU_Comm_Pack(it->first, &(ent_info[0]),3*sizeof(int));
        PCU_COMM_PACK(it->first, e);
      }
    }
  } //  for (int dim=0; dim<2; ++dim)

  PCU_Comm_Send();
  while (PCU_Comm_Listen())
  {
    while (! PCU_Comm_Unpacked())
    {
      int sender_ent_info[3];
      pMeshEnt sender;
      PCU_Comm_Unpack(&(sender_ent_info[0]), 3*sizeof(int));
      PCU_COMM_UNPACK(sender);
      e = partbdry_entities[sender_ent_info[0]][sender_ent_info[1]];
      set_remote(m, e, sender_ent_info[2], sender);
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())
}


// *********************************************************
void send_entities(pMesh m, int dim, pTag partbdry_id_tag)
// *********************************************************
{
  int own_partid, num_ent = m->count(dim);
  if (dim>2 || !num_ent) return;

  double* v_coords;
  double* v_params;
  int* e_geom_type = new int[num_ent];
  int* e_geom_tag = new int[num_ent];
  int* f_down_num; // for faces
  int* e_own_partid; // for vertices and edges
  int* e_global_id; // for vertices and edges
  int* e_rmt_num;  // for vertices and edges
  
  pMeshEnt e;
  pGeomEnt geom_ent;

  if (dim==0)
  {
    v_coords = new double[num_ent*3];  // for vertices
    v_params = new double[num_ent*3];
  }
  if (dim==2)
    f_down_num = new int[num_ent]; // for faces
  else // for vertices and edges
  {
    e_own_partid = new int[num_ent]; 
    e_global_id = new int[num_ent]; 
    e_rmt_num = new int[num_ent];  
  }

  int num_down=0, num_remote=0;
  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
  {
    switch (m->getType(e))
    {
      case apf::Mesh::TRIANGLE: num_down+=3; break;
      case apf::Mesh::QUAD: num_down+=4; break;
      case apf::Mesh::EDGE: num_down+=2; break;
      default: break;
    }

    if (m->isShared(e))
    {
      pCopies remotes;
      m->getRemotes(e, remotes);
      num_remote+=remotes.size();
    }
  }
  m->end(it);

  std::vector<int> e_down_lid;
  int* e_rmt_pid = new int[num_remote];

  int global_id, rinfo_pos=0, index=0, num_down_ent;

  it = m->begin(dim);
  apf::Vector3 coord;
  apf::Vector3 param;
  apf::Downward down_ent;
  while ((e = m->iterate(it)))
  {
    switch (getDimension(m, e))
    {
      case 0: m->getPoint(e, 0, coord);
              m->getParam(e, param);
              for (int i=0; i<3; ++i)
              {
                v_coords[index*3+i] = coord[i];
                v_params[index*3+i] = param[i];
              }
              break;
       default: { 
                  num_down_ent =  m->getDownward(e, dim-1, down_ent); 
                  if (dim==2) 
                    f_down_num[index] = num_down_ent;
                  for (int i=0; i<num_down_ent; ++i)
                    e_down_lid.push_back(getMdsIndex(m, down_ent[i]));
                }
    } // switch

    geom_ent = (gmi_ent*)(m->toModel(e));
    e_geom_type[index] = gmi_dim(pumi::instance()->model, geom_ent);
    e_geom_tag[index] = gmi_tag(pumi::instance()->model, geom_ent);

    if (dim!=2) // for vertices and edges
    {
      own_partid=m->getOwner(e);
      e_own_partid[index] = own_partid;
      pCopies remotes;
      m->getRemotes(e, remotes);
      if (m->hasTag(e,partbdry_id_tag)) // part-bdry entity
      {
        m->getIntTag(e, partbdry_id_tag, &global_id);
        e_global_id[index] = global_id; 
        e_rmt_num[index] = remotes.size();
        APF_ITERATE(pCopies, remotes, it)
        {
          e_rmt_pid[rinfo_pos]=it->first;
          ++rinfo_pos;
        }
      }
      else
      {
        e_global_id[index] = -1; 
        e_rmt_num[index] = 0;
      }
    }
    ++index;
  }
  m->end(it);

  int proc=PCU_Comm_Self()+pumi::instance()->plane_size;
  while (proc<PCU_Comm_Peers())
  {  
    PCU_Comm_Pack(proc, &num_ent, sizeof(int));
    switch (dim)
    { 
      case 0: PCU_Comm_Pack(proc, &(v_coords[0]),num_ent*3*sizeof(double));
              PCU_Comm_Pack(proc, &(v_params[0]),num_ent*3*sizeof(double));
              break;
      default: if (dim==2)
               {
                 PCU_Comm_Pack(proc, &num_down, sizeof(int));
                 PCU_Comm_Pack(proc, &(f_down_num[0]),num_ent*sizeof(int));
               }
               PCU_Comm_Pack(proc, &e_down_lid.at(0), num_down*sizeof(int));
    }
    PCU_Comm_Pack(proc, &(e_geom_type[0]),num_ent*sizeof(int));
    PCU_Comm_Pack(proc, &(e_geom_tag[0]),num_ent*sizeof(int));
    if (dim<2)  
    {      
      PCU_Comm_Pack(proc, &(e_own_partid[0]), num_ent*sizeof(int));
      PCU_Comm_Pack(proc, &(e_global_id[0]),num_ent*sizeof(int));
      PCU_Comm_Pack(proc, &(e_rmt_num[0]),num_ent*sizeof(int));
      PCU_Comm_Pack(proc, &num_remote, sizeof(int));
      if (num_remote)
        PCU_Comm_Pack(proc, &(e_rmt_pid[0]), num_remote*sizeof(int));
    }
    proc+=pumi::instance()->plane_size;
  }

  if (dim==0)
  {
    delete [] v_coords;
    delete [] v_params;
  }
  delete [] e_geom_type;
  delete [] e_geom_tag;

  if (dim==2)
    delete [] f_down_num;
  else 
  {    
    delete [] e_own_partid;
    delete [] e_global_id;
    delete [] e_rmt_num;
    if (num_remote) delete [] e_rmt_pid;
  }
}


// *********************************************************
void receive_vertices(pMesh m, pTag partbdry_id_tag, 
                           std::map<int, pMeshEnt>* partbdry_entities)
// *********************************************************
{
  int remote_rank_offset = pumi::instance()->local_planeid*pumi::instance()->plane_size;

  int myrank = PCU_Comm_Self();
  int num_ent, num_remote, own_partid;
  pMeshEnt new_ent;
  pGeomEnt geom_ent;
  while (PCU_Comm_Listen())
  {
    while (!PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      double* v_coords = new double[num_ent*3];
      double* v_params = new double[num_ent*3];
      int* e_geom_type = new int[num_ent];
      int* e_geom_tag = new int[num_ent];
      int* e_own_partid = new int[num_ent];
      int* e_global_id = new int[num_ent];
      int* e_rmt_num = new int[num_ent];

      PCU_Comm_Unpack(&(v_coords[0]), num_ent*3*sizeof(double));
      PCU_Comm_Unpack(&(v_params[0]), num_ent*3*sizeof(double));
      PCU_Comm_Unpack(&(e_geom_type[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_tag[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_own_partid[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_global_id[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_rmt_num[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&num_remote, sizeof(int));
      int* e_rmt_pid;
      if (num_remote) 
      {
        e_rmt_pid = new int[num_remote];
        PCU_Comm_Unpack(&(e_rmt_pid[0]), num_remote*sizeof(int));
      }
      // create vertex
      int rinfo_pos=0;
      apf::Vector3 coord;
      apf::Vector3 param;     
      for (int index=0; index<num_ent; ++index)
      {
        geom_ent = gmi_find(pumi::instance()->model, e_geom_type[index], e_geom_tag[index]);
        for (int k=0; k<3; ++k)
        {
          coord[k] = v_coords[index*3+k];
          param[k] = v_params[index*3+k];
        }
        new_ent = m->createVertex((apf::ModelEntity*)geom_ent, coord, param);
        if (e_global_id[index]!=-1)
        {
          partbdry_entities[0][e_global_id[index]] = new_ent;
          m->setIntTag(new_ent, partbdry_id_tag, &(e_global_id[index]));
        }

        for (int i=0; i<e_rmt_num[index]; ++i)
        {
          m->addRemote(new_ent, remote_rank_offset+e_rmt_pid[rinfo_pos], NULL);
          ++rinfo_pos;
        }
      } // for index
      delete [] v_coords;
      delete [] v_params;
      delete [] e_geom_type;
      delete [] e_geom_tag;
      delete [] e_own_partid;
      delete [] e_global_id;
      delete [] e_rmt_num;
      if (num_remote)
        delete [] e_rmt_pid;
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())
}

// *********************************************************
void receive_edges(pMesh m, pTag partbdry_id_tag, std::map<int, pMeshEnt>* partbdry_entities)
// *********************************************************
{
  int myrank=PCU_Comm_Self();
  int remote_rank_offset = pumi::instance()->local_planeid*pumi::instance()->plane_size;

  int num_ent, num_remote, own_partid;
  pMeshEnt new_ent;
  pGeomEnt geom_ent;
  apf::Downward down_ent; 
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      int* e_down_lid = new int[num_ent*2];
      int* e_geom_type = new int[num_ent];
      int* e_geom_tag = new int[num_ent];
      int* e_own_partid = new int[num_ent];
      int* e_global_id = new int[num_ent];
      int* e_rmt_num = new int[num_ent];
      PCU_Comm_Unpack(&(e_down_lid[0]), num_ent*2*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_type[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_tag[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_own_partid[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_global_id[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_rmt_num[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&num_remote, sizeof(int));
      int* e_rmt_pid = new int[num_remote];
      if (num_remote) PCU_Comm_Unpack(&(e_rmt_pid[0]), num_remote*sizeof(int));

      // create edge
      int rinfo_pos=0;
      for (int index=0; index<num_ent; ++index)
      {
        geom_ent = gmi_find(pumi::instance()->model, e_geom_type[index], e_geom_tag[index]);
        down_ent[0] =  getMdsEntity(m, 0, e_down_lid[index*2]);
        down_ent[1] =  getMdsEntity(m, 0, e_down_lid[index*2+1]);
        new_ent = m->createEntity(apf::Mesh::EDGE, (apf::ModelEntity*)geom_ent, down_ent);
        if (e_global_id[index]!=-1)
        {
          partbdry_entities[1][e_global_id[index]] = new_ent;
          m->setIntTag(new_ent, partbdry_id_tag, &(e_global_id[index]));
        }
        for (int i=0; i<e_rmt_num[index]; ++i)
        {
          m->addRemote(new_ent, remote_rank_offset+e_rmt_pid[rinfo_pos], NULL);
          ++rinfo_pos;
        }
      } // for index      
      delete [] e_down_lid;
      delete [] e_geom_type;
      delete [] e_geom_tag;
      delete [] e_own_partid;
      delete [] e_global_id;
      delete [] e_rmt_num;
      delete [] e_rmt_pid;
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())

}

// *********************************************************
void receive_faces(pMesh m)
// *********************************************************
{
  int num_ent, num_down;
  pMeshEnt new_ent;
  pGeomEnt geom_ent;
  apf::Downward down_ent;
  int myrank = PCU_Comm_Self();
  while (PCU_Comm_Listen())
  {
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_Comm_Unpack(&num_ent, sizeof(int));
      PCU_Comm_Unpack(&num_down, sizeof(int));

      int* f_down_num = new int[num_ent];
      int* e_down_lid = new int[num_down];
      int* e_geom_type = new int[num_ent];
      int* e_geom_tag = new int[num_ent];
      PCU_Comm_Unpack(&(f_down_num[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_down_lid[0]), num_down*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_type[0]), num_ent*sizeof(int));
      PCU_Comm_Unpack(&(e_geom_tag[0]), num_ent*sizeof(int));
      // create face
      int dinfo_pos=0;
      for (int index=0; index<num_ent; ++index)
      {
        geom_ent =gmi_find(pumi::instance()->model, e_geom_type[index], e_geom_tag[index]);
        num_down = f_down_num[index];
        for (int i=0; i<num_down; ++i)
        {
          down_ent[i] = getMdsEntity(m, 1, e_down_lid[dinfo_pos]);
          ++dinfo_pos;
        }
        if (num_down==3)
          new_ent= m->createEntity(apf::Mesh::TRIANGLE, (apf::ModelEntity*)geom_ent, down_ent);
        else
          new_ent= m->createEntity(apf::Mesh::QUAD, (apf::ModelEntity*)geom_ent, down_ent);
      } // for index
      delete [] f_down_num;
      delete [] e_down_lid;
      delete [] e_geom_type;
      delete [] e_geom_tag;
    } // while ( ! PCU_Comm_Unpacked())
  } // while (PCU_Comm_Listen())
}

// *********************************************************
void copy_mesh_2_slave(pMesh m)
// *********************************************************
{
  int self=PCU_Comm_Self();

  // assign uniq id to part bdry entities
  pTag partbdry_id_tag = m->createIntTag("pbdry_globid", 1);

  for (int dim=0; dim<2; ++dim)
    assign_uniq_partbdry_id(m, dim, partbdry_id_tag);

  // copy 2D mesh in process group 0 to other process groups
  std::map<int, pMeshEnt> partbdry_entities[2];

  for (int dim=0; dim<3; ++dim)
  {
    PCU_Comm_Begin();
    send_entities(m, dim, partbdry_id_tag);
    PCU_Comm_Send();
    switch (dim)
    {
      case 0: receive_vertices(m, partbdry_id_tag, partbdry_entities); break;
      case 1: receive_edges(m, partbdry_id_tag, partbdry_entities); break;
      case 2: receive_faces(m); break;
      default: break;
    }   
  }

  stitch_link(m, partbdry_id_tag, partbdry_entities);

  pMeshEnt e;
  for (int dim=0; dim<2; ++dim)
  {
    apf::MeshIterator* ent_it = m->begin(dim);
    while ((e = m->iterate(ent_it)))
    {
      if (!m->isShared(e)) continue;
      pCopies remotes;
      apf::Parts parts;
      m->getRemotes(e, remotes);
      APF_ITERATE(pCopies, remotes, it)
        parts.insert(it->first);
      parts.insert(self);
      m->setResidence (e, parts); // set pclassification
      m->removeTag(e, partbdry_id_tag);
    }
    m->end(ent_it);
  }
  m->destroyTag(partbdry_id_tag);
  
  // update partition model
  m->acceptChanges();
}

pMesh pumi_mesh_create(pGeom g, const char* filename, int num_in_part, 
                       int num_plane, const char* mesh_type)
{
  if (strcmp(mesh_type,"mds"))
  {
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid mesh type "<<mesh_type<<"\n";
    return NULL;
  }
  // set proc_group_id
  int self = PCU_Comm_Self();
  pumi::instance()->plane_size = PCU_Comm_Peers()/num_plane;
  pumi::instance()->local_planeid = self/pumi::instance()->plane_size; // divide
  pumi::instance()->prev_plane_partid = (PCU_Comm_Self()-pumi::instance()->plane_size)%PCU_Comm_Peers();
  if (pumi::instance()->prev_plane_partid<0)
    pumi::instance()->prev_plane_partid = pumi::instance()->prev_plane_partid+PCU_Comm_Peers();
  pumi::instance()->next_plane_partid = (self+pumi::instance()->plane_size)%PCU_Comm_Peers();

  if (num_in_part==1) // do static partitioning
  {
    assert(num_plane==1);
    MPI_Comm prevComm = PCU_Get_Comm();
    int num_target_part = PCU_Comm_Peers()/num_in_part;
    bool isMaster = ((PCU_Comm_Self() % num_target_part) == 0);
    pMesh m = 0;
    apf::Migration* plan = 0;   
    split_comm(num_target_part);
    if (isMaster) {
      m = apf::loadMdsMesh(g, filename);
      plan = getPlan(m, num_target_part);
    }
    merge_comm(prevComm);
    pumi::instance()->mesh = apf::repeatMdsMesh(m, g, plan, num_target_part);
  }
  else
  {
    if (num_plane>1)
    {
      int groupRank = self%pumi::instance()->plane_size; // modulo
      MPI_Comm groupComm;
      MPI_Comm_split(PCU_Get_Comm(), pumi::instance()->local_planeid, groupRank, &groupComm);
      PCU_Switch_Comm(groupComm);
    }

    if (!pumi::instance()->local_planeid) // master plane
        pumi::instance()->mesh = apf::loadMdsMesh(g, filename);
    else // slave plane
      pumi::instance()->mesh = apf::makeEmptyMdsMesh(g, 2, false);

    if (num_plane>1)
    {
      merge_comm(MPI_COMM_WORLD);
      copy_mesh_2_slave(pumi::instance()->mesh);
    }
  }
  generate_node_global_numbering(pumi::instance()->mesh);
  return pumi::instance()->mesh;
}

int pumi_mesh_getdim(pMesh m)
{
  return m->getDimension();
}

int pumi_mesh_getnument(pMesh m, int dim)
{ return m->count(dim); }

#include <parma.h>
void pumi_mesh_print (pMesh m)
{
  if (!PCU_Comm_Self()) std::cout<<"\n=== mesh size and tag info === \nglobal ";
  printStats(m);

  int* local_entity_count = new int[4*PCU_Comm_Peers()];
  for (int i=0; i<4*PCU_Comm_Peers();++i)
    local_entity_count[i]=0;
  for (int d=0; d<4;++d)
    local_entity_count[4*pumi_rank()+d] = m->count(d);
  
  int* global_entity_count = new int[4*PCU_Comm_Peers()]; 

  MPI_Allreduce(local_entity_count, global_entity_count, 4*PCU_Comm_Peers(), MPI_INT, MPI_SUM, PCU_Get_Comm());
  if (pumi_rank()) return;

  for (int p=0; p<PCU_Comm_Peers(); ++p)
    std::cout<<"(p"<<p<<") # local ent: v "<<global_entity_count[p*4]
        <<", e "<<global_entity_count[p*4+1]
        <<", f "<<global_entity_count[p*4+2]
        <<", r "<<global_entity_count[p*4+3]<<"\n";

  delete [] local_entity_count;
  delete [] global_entity_count;

  apf::DynamicArray<pTag> tags;
  m->getTags(tags);
  int n = tags.getSize();
  for (int i = 0; i < n; ++i) 
    std::cout<<"tag "<<i<<": \""<< m->getTagName(tags[i])<<"\", type "<< m->getTagType(tags[i])<<", size "<< m->getTagSize(tags[i])<<"\n";
}

void pumi_mesh_write (pMesh m, const char* filename, const char* mesh_type)
{
  if (!strcmp(mesh_type,"mds"))
    m->writeNative(filename);
  else if (!strcmp(mesh_type,"vtk"))
    apf::writeVtkFiles(filename, m);
  else
    if (!PCU_Comm_Self()) std::cout<<"[PUMI ERROR] "<<__func__<<" failed: invalid mesh type "<<mesh_type<<"\n";
}

void pumi_mesh_delete(pMesh m)
{
  destroy_node_global_numbering(m);
  m->destroyNative();
  apf::destroyMesh(m);
}

void pumi_mesh_verify(pMesh m)
{
  apf::verify(m);
}
