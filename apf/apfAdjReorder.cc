/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfField.h"
#include "apfMesh.h"
#include "apfNumbering.h"
#include "apfShape.h"

#include <list>

namespace apf {

  int NaiveOrder(Numbering * num)
  {
    Field * field = getField(num);
    Mesh * mesh = getMesh(field);
    FieldShape * shape = getShape(field);
    
    int components = countComponents(field);
    int dim = mesh->getDimension();
    
    long dof = 0;

    MeshIterator * it;
    MeshEntity * ent = NULL;
    for(int ii = 0; ii < dim; ii++)
    {
      if(shape->hasNodesIn(ii))
      {
        it = mesh->begin(ii);
        while((ent = mesh->iterate(it)))
        {
          if(mesh->isOwned(ent))
          {
            int type = mesh->getType(ent);
            int type_nodes = shape->countNodesOn(type);

            for(int jj = 0; jj < type_nodes; jj++)
              for(int kk = 0; kk < components; kk++)
              {
                if(!isFixed(num,ent,jj,kk))
                {
                  number(num,ent,jj,kk,dof);
                  dof++;
                }
              }
          }
        }
        mesh->end(it);
      }
    }
    return dof;
  }


int AdjReorder(Numbering * num)
{
  Field * field = getField(num);
  Mesh * mesh = getMesh(field);
  FieldShape * shape = getShape(field);  
  int dim = mesh->getDimension();
  int components = countComponents(field);
  int dofs = 0;

  MeshIterator * it;
  MeshEntity * ent = NULL;

  // count dofs on all locally owned mesh entities
  for(int ii = 0; ii < dim; ii++)
  {
    if(shape->hasNodesIn(ii))
    {
      it = mesh->begin(ii);
      while((ent = mesh->iterate(it)))
      {
        if(mesh->isOwned(ent))
          dofs += shape->countNodesOn(mesh->getType(ent)) * components;
      }
      mesh->end(it);
    }
  }
  
  // deduct any fixed entities from the locally-owned total
  int num_fixed = countFixed(num);
  dofs -= num_fixed;

// above here works

  std::list<MeshEntity*> ent_queue;
  std::set<MeshEntity*> queue_set;

  // find the first locally-owned vertex
  it = mesh->begin(0);
  do
  {
    ent = mesh->iterate(it);
  } while(!mesh->isOwned(ent));
  ent_queue.push_back(ent);
  queue_set.insert(ent);
  mesh->end(it);

// above here works

  // on the off chance there are multiple nodes per vertex...
  int vert_nodes = shape->countNodesOn(Mesh::VERTEX);
  long dof_index = dofs-1;

  while(dof_index >= 0)
  {
    if(ent_queue.empty())
    {
      NaiveOrder(num);
      break;
    }

    ent = ent_queue.front(); ent_queue.pop_front();
    queue_set.erase(ent);

    int type = mesh->getType(ent);
    int type_nodes = shape->countNodesOn(type);

    for(int ii = 0; ii < type_nodes; ii++)
      for(int jj = 0; jj < components; jj++)
      {
        if(!isFixed(num,ent,ii,jj) && 
	   !isNumbered(num,ent,ii,jj) )
        {
          number(num,ent,ii,jj,dof_index);
	  //std::cout << dof_index << std::endl;
          dof_index--;
        }
      }

    if(type != Mesh::VERTEX)
      continue;
    else // type is Mesh::VERTEX
    {
      for(int ii = 0; ii < mesh->countUpward(ent); ii++)
      {
	MeshEntity * e = mesh->getUpward(ent,ii);

	Adjacent elmnts;
	mesh->getAdjacent(e,dim,elmnts);
	APF_ITERATE(Adjacent,elmnts,eit)
	{
	  Adjacent verts;
	  mesh->getAdjacent(*eit,0,verts);
	  APF_ITERATE(Adjacent,verts,vit)
	  {
	    if( !isNumbered(num,*vit,0,0) &&
		queue_set.find(*vit) == queue_set.end() &&
		mesh->isOwned(*vit) )
	    {
	      ent_queue.push_back(*eit);
	      queue_set.insert(*eit);
	    }
	  }
	}

	MeshEntity* ov = getEdgeVertOppositeVert(mesh,e,ent);
	bool owned = mesh->isOwned(ov);
	if(owned)
	{
	  bool queued = queue_set.find(ov) != queue_set.end();
	  bool numbered = false;
	  for(int jj = 0; jj < vert_nodes; jj++)
	    for(int kk = 0; kk < components; kk++)
	      numbered = numbered == true ? true : isNumbered(num,ov,jj,kk);

	  if(!numbered)
	  {
	    if(queued)
	    {
	      for(int jj = 0; jj < vert_nodes; jj++)
		for(int kk = 0; kk < components; kk++)
		{
		  if(!isFixed(num,ov,jj,kk))
		  {
		    number(num,ov,jj,kk,dof_index);
		    //std::cout << dof_index << std::endl;
		    dof_index--;
		  }
		}
	      ent_queue.remove(ov);
	      queue_set.erase(ov);
	    }
	    else
	    {
	      ent_queue.push_back(ov);
	      queue_set.insert(ov);
	    }
	  }
	}
      }
    }
  }
      
  return dofs;
}

void reorderConnected(
    Mesh* mesh,
    MeshTag* numbering,
    int& node_label,
    int& element_label)
{
  int dimension = mesh->getDimension();
  MeshEntity* startVertex = 0;
  MeshEntity* v;
  int bestModelDimension = dimension+1;
  MeshIterator* vertices = mesh->begin(0);
  while ((v = mesh->iterate(vertices)))
  {
    if (mesh->hasTag(v,numbering))
      continue;
    int modelDimension = mesh->getModelType(mesh->toModel(v));
    if (modelDimension < bestModelDimension)
    {
      startVertex = v;
      bestModelDimension = modelDimension;
    }
  }
  mesh->end(vertices);
  std::list<MeshEntity*> queue;
  std::set<MeshEntity*> vset;
  queue.push_back(startVertex);
  vset.insert(startVertex);
  while (( ! queue.empty())&&(node_label >= 0))
  {
    v = queue.front(); queue.pop_front();
    vset.erase(v);
    if ( ! mesh->hasTag(v,numbering))
    {
      mesh->setIntTag(v,numbering,&node_label);
      --node_label;
    }
    for (int i=0; i < mesh->countUpward(v); ++i)
    {
      MeshEntity* e = mesh->getUpward(v,i);
      MeshEntity* ov = getEdgeVertOppositeVert(mesh,e,v);
      if (( ! mesh->hasTag(ov,numbering))
        &&( vset.find(ov) == vset.end() ))
      {
        queue.push_back(ov);
        vset.insert(ov);
      }
      Adjacent elements;
      mesh->getAdjacent(e,dimension,elements);
      APF_ITERATE(Adjacent,elements,eit)
      {
        if ( ! mesh->hasTag(*eit,numbering))
        {
          mesh->setIntTag(*eit,numbering,&element_label);
          --element_label;
        }
      }
    }
  }
}

MeshTag* reorder(Mesh* mesh, const char* name)
{
  int dimension = mesh->getDimension();
  int node_label = mesh->count(0)-1;
  int element_label = mesh->count(dimension)-1;
  MeshTag* numbering = mesh->createIntTag(name,1);
  while (node_label >= 0)
    reorderConnected(mesh,numbering,node_label,element_label);
  assert(node_label==-1);
  assert(element_label==-1);
  return numbering;
}

void SetNumberingOffset(Numbering * num, int off)
{
  Field * field = getField(num);
  Mesh * mesh = getMesh(field);
  FieldShape * shape = getShape(field);

  int components = countComponents(field);
  int dim = mesh->getDimension();

  /* iterate over all nodes in the mesh,
     get their current numbering,
     add the offset and set the new numbering */
  MeshIterator * iter;
  MeshEntity * e;
  for(int ii = 0; ii < dim; ii++)
  {
    if(!shape->hasNodesIn(ii))
      break;

    iter = mesh->begin(ii);
    while((e = mesh->iterate(iter)))
    {
      if(mesh->isOwned(e))
      {
	int type = mesh->getType(e);
	int type_nodes = shape->countNodesOn(type);
	
	for(int ii = 0; ii < type_nodes; ii++)
	  for(int jj = 0; jj < components; jj++)
	  {
	    if(isNumbered(num,e,ii,jj))
	    {
	      int current = getNumber(num,e,ii,jj);
	      number(num,e,ii,jj,current+off);
	      //std::cout << current+off << std::endl;
	    }
	  }
      }
    }
    mesh->end(iter);
  }
}

}

