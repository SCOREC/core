/****************************************************************************** 

  (c) 2004-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef MPART_ENTITY_CONTAINER_H
#define MPART_ENTITY_CONTAINER_H

#include "GenTag.h"
#include "pumi_list.h"
#include "gmi.h"

class gEntity : public Taggable, public ListMember
{
public:
  gEntity(gmi_ent*);
  ~gEntity();
  gmi_ent* getGmi() {return e;}  
private:
  gmi_ent* e;
};

class mPartEntityContainer
{
  public :
    typedef List* CONTAINER;
    typedef ListIterator<gEntity> iter;
  private:
    enum { _DIMS_ = 4 };
    List gEntities[_DIMS_];
  public:
    mPartEntityContainer();
    virtual ~mPartEntityContainer();
    iter begin(int what);
    iter end(int what);
    void add(int d, gEntity* e);
    void del(int d, gEntity* e);
    int  size(int what) const;
};

#endif
