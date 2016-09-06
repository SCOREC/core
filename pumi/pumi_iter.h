/****************************************************************************** 

  (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef _PUMI_ITERATOR_H_
#define _PUMI_ITERATOR_H_
#include "pumi.h"
#include "FMDBfwd.h"
#include "mEntity.h"
#include "GenIterator.h"
#include "mPart.h"
#include "mEntitySet.h"
#include "pumi_mesh.h"

typedef mPartEntityContainer::iter part_iter;

typedef EntitySetUnordered<mEntity>::iter entSetU_iter;
typedef EntitySetOrdered<mEntity>::iter entSetO_iter;

struct part_pid_struct
{
   void* part;
   int pid;
};

struct part_gent_struct
{
  void* part;
  void* gent;
};


inline void processingEntitySetOFilter(entSetO_iter& it_begin, entSetO_iter& it_end, void* ptr, int type, int topo)
{
  if(it_begin==it_end)
    return; 
  entSetO_iter it;
  
  switch (topo)
  {
    case PUMI_POINT:
    case PUMI_LINE:
      for (it= it_begin; it!=it_end;++it)
	 if ((*it)!=NULL && (*it)->getTopo()==topo) {
	    it_begin = it; 
	    return; 
	 }  
	 break;
    case PUMI_ALLTOPO:
      if (type == 4) // ALL_TYPES
      {
	 return; 
      }
      else // SPECIFIC TYPE
      {
        for (it=it_begin; it!=it_end;++it)
          if((*it)->getLevel()==type) {
            it_begin = it; 
	    return;
	  }
      }
      break;
    default: if (type == PUMI_ALLTYPE)
           {
             for (it=it_begin; it!=it_end;++it)
               if ((*it)->getTopo()==topo) {
                 it_begin = it; 
		 return; 
	       }
           }
           else // SPECIFIC TYPE & SPECIFIC TOPO
           {
             for (it=it_begin; it!=it_end;++it)
               if ((*it)->getLevel()==type && (*it)->getTopo()==topo) {
                 it_begin = it; 
		 return; 
	       }
           }
           break;
  } // end of switch
  it_begin=it_end; 
}

inline void processingEntitySetUFilter(entSetU_iter& it_begin, entSetU_iter& it_end, void* ptr, int type, int topo)
{
  if(it_begin==it_end)
    return; 
  entSetU_iter it;
  switch (topo)
  {
    case PUMI_POINT:
    case PUMI_LINE: for (it= it_begin; it!=it_end;++it)
                                          if ((*it)->getTopo()==topo) {
                                             it_begin = it; 
					     return;
					  }
                                      break;
    case PUMI_ALLTOPO:
      if (type == PUMI_ALLTYPE)
      {
        return; 
      }
      else // SPECIFIC TYPE
      {
        for (it=it_begin; it!=it_end;++it)
          if((*it)->getLevel()==type) {
            it_begin = it; 
	    return; 
	  }
      }
      break;
    default: if (type == PUMI_ALLTYPE)
           {
             for (it=it_begin; it!=it_end;++it)
               if ((*it)->getTopo()==topo) {
                 it_begin = it; 
		 return;
	       }
           }
           else // SPECIFIC TYPE & SPECIFIC TYPE
           {
             for (it=it_begin; it!=it_end;++it)
               if ((*it)->getLevel()==type && (*it)->getTopo()==topo) {
                 it_begin = it; 
		 return; 
	       }
           }
           break;
  } // end of switch
  it_begin=it_end; 
}
 
void SingleTypeFilter(part_iter& it_begin, part_iter& it_end);
void SingleTypeTopoFilter(part_iter& it_begin, part_iter& it_end, int topo);
void AllTypeFilter(part_iter& it_begin, part_iter& it_end, pPart part);

inline void processingMeshFilter(part_iter& it_begin, part_iter& it_end, void* ptr, int type, int topo)
{
  if (it_begin==it_end) return;

  mPart* part = (mPart*)ptr;    
  switch (type)
  {     
    case PUMI_VERTEX:
    case PUMI_EDGE: SingleTypeFilter(it_begin, it_end); break;
    case PUMI_FACE: if (topo==PUMI_ALLTOPO)
                      SingleTypeFilter(it_begin, it_end);
                    else
                      SingleTypeTopoFilter(it_begin, it_end, topo);
                    break;
    case PUMI_REGION: if (topo==PUMI_ALLTOPO)
                        SingleTypeFilter(it_begin, it_end);
                      else
                        SingleTypeTopoFilter(it_begin, it_end, topo);  
                      break;
    case PUMI_ALLTYPE: AllTypeFilter(it_begin, it_end, part);
                    break;
    default: break; // all type & specific topo pair is impossible -- see "api/pumi_iter.cc" 285L
  } 
}

#ifdef PUMI_PARALLEL
void PartBdrySingleTypeFilter(part_iter& it_begin, part_iter& it_end, int target_pid);
void PartBdrySingleTypeTopoFilter(part_iter& it_begin, part_iter& it_end, int target_pid, int topo);
void PartBdryAllTypeFilter(part_iter& it_begin, part_iter& it_end, pPart part, int target_pid);
#endif

inline void processingPartBdryEntFilter(part_iter& it_begin, part_iter& it_end, void* ptr, int type, int topo)
{
#ifdef PUMI_PARALLEL
  if (it_begin==it_end || type==PUMI_REGION || PUMI_Topo_Type(topo)==3) 
  {  it_begin = it_end; return; }

  part_pid_struct* part_pid = (part_pid_struct*)ptr; 
  pPart part = (pPart)(part_pid->part); 
  int target_pid = part_pid->pid; 

  switch (type)
  {
    case PUMI_VERTEX: 
    case PUMI_EDGE: PartBdrySingleTypeFilter(it_begin, it_end, target_pid); break;
    case PUMI_FACE: if (topo==PUMI_ALLTOPO) 
                      PartBdrySingleTypeFilter(it_begin, it_end, target_pid); 
                    else
                      PartBdrySingleTypeTopoFilter(it_begin, it_end, target_pid, topo); 
                    break;
    case PUMI_ALLTYPE: PartBdryAllTypeFilter(it_begin, it_end, part, target_pid); 
                    break;
    default: break; // all type & specific topo pair is impossible -- see "api/pumi_iter.cc" 285L
  }
#endif
}

void RevClasSingleTypeFilter(part_iter& it_begin, part_iter& it_end, pGeomEnt gent);
void RevClasAllTypeFilter(part_iter& it_begin, part_iter& it_end, pPart part, pGeomEnt gent);

inline void processingGeomFilter(part_iter& it_begin, part_iter& it_end, void* ptr, int type, int topo)
{
  if (it_begin==it_end)
    return; 

  part_gent_struct* part_gent = (part_gent_struct*)ptr; 
  pPart part = (pPart)(part_gent->part);
  pGeomEnt gent = (pGeomEnt)(part_gent->gent);

  if (type==PUMI_ALLTYPE)
    RevClasAllTypeFilter(it_begin, it_end, part, gent);
  else
    RevClasSingleTypeFilter(it_begin, it_end, gent);
}

inline void processingNonFilter(part_iter& it_begin, part_iter& it_end, void* ptr, int type, int topo)
{}

#endif
