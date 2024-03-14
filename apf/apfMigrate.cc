/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfMesh2.h"
#include "apfCavityOp.h"
#include "apf.h"
#include <pcu_util.h>
#include <lionPrint.h>
#include <cstdlib>

namespace apf {

typedef std::vector<MeshEntity*> EntityVector;

/* Starting from the elements in the plan,
   constructs their closure (including all
   remote copies of the closure). */
/* if there is matching, bring all the matches
   of an affected entity into the affected set as well. */
static void getAffected(
    Mesh2* m,
    Migration* plan,
    EntityVector affected[4])
{
  int maxDimension = m->getDimension();
  int self = m->getPCU()->Self();
  affected[maxDimension].reserve(plan->count());
  for (int i=0; i < plan->count(); ++i)
  {
    MeshEntity* e = plan->get(i);
    if (plan->sending(e) != self) {
      PCU_ALWAYS_ASSERT(apf::getDimension(m, e) == m->getDimension());
      affected[maxDimension].push_back(e);
    }
  }
  int dummy;
  MeshTag* tag = m->createIntTag("apf_migrate_affected",1);
  for (int dimension=maxDimension-1; dimension >= 0; --dimension)
  {
    int upDimension = dimension + 1;
    m->getPCU()->Begin();
    APF_ITERATE(EntityVector,affected[upDimension],it)
    {
      MeshEntity* up = *it;
      Downward adjacent;
      int na = m->getDownward(up,dimension,adjacent);
      for (int i=0; i < na; ++i)
      {
        if ( ! m->hasTag(adjacent[i],tag))
        {
          m->setIntTag(adjacent[i],tag,&dummy);
          affected[dimension].push_back(adjacent[i]);
        }
        Copies remotes;
        m->getRemotes(adjacent[i],remotes);
        APF_ITERATE(Copies,remotes,rit)
          m->getPCU()->Pack(rit->first,rit->second);
        if (m->hasMatching())
        {
          Matches matches;
          m->getMatches(adjacent[i],matches);
          for (size_t j=0; j < matches.getSize(); ++j)
            m->getPCU()->Pack(matches[j].peer,matches[j].entity);
        }
      }//downward adjacent loop
    }//upward affected loop
    m->getPCU()->Send();
    while (m->getPCU()->Receive())
    {
      MeshEntity* entity;
      m->getPCU()->Unpack(entity);
      if ( ! m->hasTag(entity,tag))
      {
        m->setIntTag(entity,tag,&dummy);
        affected[dimension].push_back(entity);
      }
    }
    APF_ITERATE(EntityVector,affected[dimension],it)
      m->removeTag(*it,tag);
  }//dimension loop
  m->destroyTag(tag);
}

/* gets the subset of the closure copies
   which are owned by this part */
void getSenders(
    Mesh2* m,
    EntityVector affected[4],
    EntityVector senders[4])
{
  for (int i=0; i < 4; ++i)
  {
    /* maybe overkill for pre-allocation, but
       ensures the vector will never re-allocate. */
    senders[i].reserve(affected[i].size());
    APF_ITERATE(EntityVector,affected[i],it)
      if (m->isOwned(*it))
        senders[i].push_back(*it);
  }
}

/* at this point if one matched copy is affected, all of
   its matches are in the affected set and the "senders", or
   representatives of a set of remote copies, are in
   the senders vector.
   So, if each sender just remembers the other senders
   it is matched to, that is enough for them to re-negotiate
   matches after each of them does the job of creating new
   entities and computing new remote copies */
void reduceMatchingToSenders(
    Mesh2* m,
    EntityVector senders[4])
{
  if ( ! m->hasMatching()) return;
  m->getPCU()->Begin();
  for (int d=0; d < 4; ++d)
  {
    for (size_t i=0; i < senders[d].size(); ++i)
    {
      MeshEntity* e = senders[d][i];
      Matches matches;
      m->getMatches(e,matches);
      for (size_t j=0; j < matches.getSize(); ++j)
      { /* advertise to the match that this is a sender */
        int to = matches[j].peer;
        m->getPCU()->Pack(to,matches[j].entity);
        m->getPCU()->Pack(to,e);
      }
      /* now forget all other matchings */
      m->clearMatches(e);
    }
  }
  m->getPCU()->Send();
  while (m->getPCU()->Listen())
  {
    int sender = m->getPCU()->Sender();
    while ( ! m->getPCU()->Unpacked())
    {
      MeshEntity* e;
      m->getPCU()->Unpack(e);
      MeshEntity* match;
      m->getPCU()->Unpack(match);
/* all matches of the sender receive this.
   We only care that the senders themselves receive
   messages from only other senders matched to them,
   but some non-senders do harmlessly add a few matchings */
      m->addMatch(e,sender,match);
    }
  }
/* now senders only match other senders */
}

static Parts makeResidence(int part)
{
  Parts r;
  r.insert(part);
  return r;
}

void packParts(int to, Parts& parts, pcu::PCU *PCUObj)
{
  size_t n = parts.size();
  PCUObj->Pack(to,n);
  APF_ITERATE(Parts,parts,it)
  {
    int p = *it;
    PCUObj->Pack(to,p);
  }
}

void unpackParts(Parts& parts, pcu::PCU *PCUObj)
{
  size_t n;
  PCUObj->Unpack(n);
  for (size_t i=0;i<n;++i)
  {
    int p;
    PCUObj->Unpack(p);
    parts.insert(p);
  }
}

/* for every entity in the affected closure,
   this function changes the residence to be
   the union of all upward adjacent residences
   (including those of remote copies) */
static void updateResidences(
    Mesh2* m,
    Migration* plan,
    EntityVector affected[4])
{
  int maxDimension = m->getDimension();
  for (int i=0; i < plan->count(); ++i)
  {
    MeshEntity* e = plan->get(i);
    Parts res = makeResidence(plan->sending(e));
    m->setResidence(e,res);
  }
  for (int dimension = maxDimension-1; dimension >= 0; --dimension)
  {
    m->getPCU()->Begin();
    APF_ITERATE(EntityVector,affected[dimension],it)
    {
      MeshEntity* entity = *it;
      Parts newResidence;
      Up upward;
      m->getUp(entity, upward);
      for (int ui=0; ui < upward.n; ++ui)
      {
        MeshEntity* up = upward.e[ui];
        Parts upResidence;
        m->getResidence(up,upResidence);
        unite(newResidence,upResidence);
      }
      m->setResidence(entity,newResidence);
      Copies remotes;
      m->getRemotes(entity,remotes);
      APF_ITERATE(Copies,remotes,rit)
      {
        m->getPCU()->Pack(rit->first,rit->second);
        packParts(rit->first,newResidence, m->getPCU());
      }
    }
    m->getPCU()->Send();
    while(m->getPCU()->Receive())
    {
      MeshEntity* entity;
      m->getPCU()->Unpack(entity);
      Parts current;
      m->getResidence(entity,current);
      Parts incoming;
      unpackParts(incoming, m->getPCU());
      unite(current,incoming);
      m->setResidence(entity,current);
    }
  }
}

/* given the sets of old remote copies
   and new resident parts, this function
   constructs the set
   of REMOTE parts in the new residence that
   don't have remotes yet, which is the
   set we have to send to */
void split(
    Copies& remotes,
    Parts& parts,
    Parts& newParts,
    pcu::PCU *PCUObj)
{
  APF_ITERATE(Parts,parts,it)
    if (( ! remotes.count(*it))&&(*it != PCUObj->Self()))
      newParts.insert(*it);
}

static void packCommon(
    Mesh2* m,
    int to,
    MeshEntity* e)
{
  m->getPCU()->Pack(to,e);
  ModelEntity* me = m->toModel(e);
  int modelType = m->getModelType(me);
  m->getPCU()->Pack(to,modelType);
  int modelTag = m->getModelTag(me);
  m->getPCU()->Pack(to,modelTag);
  Parts residence;
  m->getResidence(e,residence);
  packParts(to,residence, m->getPCU());
}

void unpackCommon(
    Mesh2* m,
    MeshEntity*& sender,
    ModelEntity*& c,
    Parts& residence)
{
  m->getPCU()->Unpack(sender);
  int modelType,modelTag;
  m->getPCU()->Unpack(modelType);
  m->getPCU()->Unpack(modelTag);
  c = m->findModelEntity(modelType,modelTag);
  unpackParts(residence, m->getPCU());
}

static void packVertex(
    Mesh2* m,
    int to,
    MeshEntity* e)
{
  Vector3 p;
  m->getPoint(e,0,p);
  m->getPCU()->Pack(to,p);
  m->getParam(e,p);
  m->getPCU()->Pack(to,p);
}

MeshEntity* unpackVertex(
    Mesh2* m,
    ModelEntity* c)
{
  Vector3 point;
  m->getPCU()->Unpack(point);
  Vector3 param;
  m->getPCU()->Unpack(param);
  return m->createVertex(c,point,param);
}

static void packReference(
    Mesh2* m,
    int to,
    MeshEntity* e)
{
  Copies remotes;
  m->getRemotes(e,remotes);
  Copies::iterator found = remotes.find(to);
  if (found!=remotes.end())
  {
    MeshEntity* remote = found->second;
    m->getPCU()->Pack(to,remote);
  }
  else
  {
    Copies ghosts;
    m->getGhosts(e,ghosts);
    found = ghosts.find(to);
    PCU_ALWAYS_ASSERT(found!=ghosts.end());
    MeshEntity* ghost = found->second;
    m->getPCU()->Pack(to,ghost);
  }
}

static void packDownward(Mesh2* m, int to, MeshEntity* e)
{
  Downward down;
  int d = getDimension(m, e);
  int n = m->getDownward(e,d-1,down);
  m->getPCU()->Pack(to,n);
  for (int i=0; i < n; ++i)
    packReference(m,to,down[i]);
}

static void unpackDownward(
    Downward& entities,
    pcu::PCU *PCUObj)
{
  int n;
  PCUObj->Unpack(n);
  for (int i=0; i < n; ++i)
    PCUObj->Unpack(entities[i]);
}

static void packNonVertex(
    Mesh2* m,
    int to,
    MeshEntity* e)
{
  packDownward(m,to,e);
}

MeshEntity* unpackNonVertex(
    Mesh2* m,
    int type, ModelEntity* c)
{
  Downward down;
  unpackDownward(down, m->getPCU());
  return m->createEntity(type,c,down);
}

static void packTags(
    Mesh2* m,
    int to,
    MeshEntity* e,
    DynamicArray<MeshTag*>& tags)
{
  size_t total = tags.getSize();
  size_t n = 0;
  for (size_t i=0; i < total; ++i)
    if (m->hasTag(e,tags[i]))
      ++n;
  m->getPCU()->Pack(to,n);
  for (size_t i=0; i < total; ++i)
  {
    MeshTag* tag = tags[i];
    if (m->hasTag(e,tag))
    {
      m->getPCU()->Pack(to,i);
      int type = m->getTagType(tag);
      int size = m->getTagSize(tag);
      if (type == Mesh2::DOUBLE)
      {
        DynamicArray<double> d(size);
        m->getDoubleTag(e,tag,&(d[0]));
        m->getPCU()->Pack(to,&(d[0]),size*sizeof(double));
      }
      if (type == Mesh2::INT)
      {
        DynamicArray<int> d(size);
        m->getIntTag(e,tag,&(d[0]));
        m->getPCU()->Pack(to,&(d[0]),size*sizeof(int));
      }
    }
  }
}

// seol
void packRemotes(
    Mesh2* m,
    int to,
    MeshEntity* e)
{
  Copies remotes;
  m->getRemotes(e,remotes);
  size_t n = remotes.size();
  m->getPCU()->Pack(to,n);
  APF_ITERATE(Copies,remotes,rit)
  {
    int p=rit->first;
    MeshEntity* remote = rit->second;
    m->getPCU()->Pack(to,p);
    m->getPCU()->Pack(to,remote);
  }
}

void unpackRemotes(Mesh2* m, MeshEntity* e)
{
  size_t n;
  m->getPCU()->Unpack(n);
  for (size_t i=0; i < n; ++i)
  {
    int p;
    m->getPCU()->Unpack(p);
    MeshEntity* r;
    m->getPCU()->Unpack(r);
    m->addRemote(e, p, r);
  }
}

void unpackTags(
    Mesh2* m,
    MeshEntity* e,
    DynamicArray<MeshTag*>& tags)
{
  size_t n;
  m->getPCU()->Unpack(n);
  PCU_ALWAYS_ASSERT_VERBOSE(n<=tags.size(),
      "A tag was created that does not exist on all processes.");
  for (size_t t=0; t < n; ++t)
  {
    size_t i;
    m->getPCU()->Unpack(i);
    MeshTag* tag = tags[i];
    int type = m->getTagType(tag);
    int size = m->getTagSize(tag);
    if (type == Mesh2::DOUBLE)
    {
      DynamicArray<double> d(size);
      m->getPCU()->Unpack(&(d[0]),size*sizeof(double));
      m->setDoubleTag(e,tag,&(d[0]));
    }
    if (type == Mesh2::INT)
    {
      DynamicArray<int> d(size);
      m->getPCU()->Unpack(&(d[0]),size*sizeof(int));
      m->setIntTag(e,tag,&(d[0]));
    }
  }
}

void packEntity(
    Mesh2* m,
    int to,
    MeshEntity* e,
    DynamicArray<MeshTag*>& tags,
    bool ghosting)
{
  int type = m->getType(e);
  m->getPCU()->Pack(to,type);
  packCommon(m,to,e);
  if (type == Mesh::VERTEX)
    packVertex(m,to,e);
  else
    packNonVertex(m,to,e);
  packTags(m,to,e,tags);
  if (ghosting) 
    packRemotes(m, to, e);
}

static MeshEntity* unpackEntity(
    Mesh2* m,
    DynamicArray<MeshTag*>& tags)
{
  int from = m->getPCU()->Sender();
  int type;
  m->getPCU()->Unpack(type);
  MeshEntity* sender;
  ModelEntity* c;
  Parts residence;
  unpackCommon(m,sender,c,residence);
  MeshEntity* entity;
  if (type == Mesh::VERTEX)
    entity = unpackVertex(m,c);
  else
    entity = unpackNonVertex(m,type,c);
  m->setResidence(entity,residence);
  unpackTags(m,entity,tags);
  Copies remotes;
  /* temporarily store the sender as
     the only remote copy */
  m->addRemote(entity, from, sender);
  return entity;
}

static void sendEntities(
    Mesh2* m,
    EntityVector& senders,
    DynamicArray<MeshTag*>& tags)
{
  APF_ITERATE(EntityVector,senders,it)
  {
    MeshEntity* entity = *it;
    Copies remotes;
    m->getRemotes(entity,remotes);
    Parts residence;
    m->getResidence(entity,residence);
    Parts sendTo;
    split(remotes,residence,sendTo,m->getPCU());
    APF_ITERATE(Parts,sendTo,sit)
      packEntity(m,*sit,entity,tags);
  }
}

static void receiveEntities(
    Mesh2* m,
    DynamicArray<MeshTag*>& tags,
    EntityVector& received)
{
  received.reserve(1024);
  while (m->getPCU()->Receive())
    received.push_back(unpackEntity(m,tags));
}

static void echoRemotes(
    Mesh2* m,
    EntityVector& received)
{
  APF_ITERATE(EntityVector,received,it)
  {
    MeshEntity* entity = *it;
    /* unpackEntity() stored the identity of the
       sender as the only remote copy so we could use
       it here to echo copies back to their sender */
    Copies temp;
    m->getRemotes(entity,temp);
    int from = temp.begin()->first;
    MeshEntity* sender = temp.begin()->second;
    m->getPCU()->Pack(from,sender);
    m->getPCU()->Pack(from,entity);
  }
}

static void receiveRemotes(Mesh2* m)
{
  while (m->getPCU()->Listen())
  {
    int from = m->getPCU()->Sender();
    while ( ! m->getPCU()->Unpacked())
    {
      MeshEntity* sender;
      m->getPCU()->Unpack(sender);
      MeshEntity* entity;
      m->getPCU()->Unpack(entity);
      PCU_ALWAYS_ASSERT(entity);
      m->addRemote(sender, from, entity);
    }
  }
}

/* at this stage we have all old and
   new remote copies. this function
   filters out the set of new copies,
   which is the set of remotes in
   the new residence plus this copy
   if its also in the new residence */
static void getNewCopies(
    Mesh2* m,
    MeshEntity* e,
    Copies& allRemotes,
    Copies& newCopies)
{
  Parts residence;
  m->getResidence(e,residence);
  APF_ITERATE(Copies,allRemotes,it)
    if (residence.count(it->first))
      newCopies.insert(*it);
  int rank = m->getPCU()->Self();
  if (residence.count(rank))
    newCopies[rank]=e;
}

static void packCopies(
    int to,
    Copies& copies,
    pcu::PCU *PCUObj)
{
  int n = copies.size();
  PCUObj->Pack(to,n);
  APF_ITERATE(Copies,copies,it)
  {
    PCUObj->Pack(to,it->first);
    PCUObj->Pack(to,it->second);
  }
}

static void unpackCopies(
    Copies& copies,
    pcu::PCU *PCUObj)
{
  int n;
  PCUObj->Unpack(n);
  for (int i=0; i < n; ++i)
  {
    int part;
    PCUObj->Unpack(part);
    MeshEntity* remote;
    PCUObj->Unpack(remote);
    PCU_ALWAYS_ASSERT(remote);
    copies[part]=remote;
  }
}

static void bcastRemotes(
    Mesh2* m,
    EntityVector& senders)
{
  m->getPCU()->Begin();
  int rank = m->getPCU()->Self();
  APF_ITERATE(EntityVector,senders,it)
  {
    MeshEntity* e = *it;
    Copies allRemotes;
    m->getRemotes(e,allRemotes);
    Copies newCopies;
    getNewCopies(m,e,allRemotes,newCopies);
    APF_ITERATE(Copies,allRemotes,rit)
    {
      m->getPCU()->Pack(rit->first,rit->second);
      packCopies(rit->first,newCopies, m->getPCU());
    }
    newCopies.erase(rank);
    m->setRemotes(e,newCopies);
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive())
  {
    MeshEntity* e;
    m->getPCU()->Unpack(e);
    Copies copies;
    unpackCopies(copies, m->getPCU());
    copies.erase(rank);
    m->setRemotes(e,copies);
  }
}

static void setupRemotes(
    Mesh2* m,
    EntityVector& received,
    EntityVector& senders)
{
  m->getPCU()->Begin();
  echoRemotes(m,received);
  m->getPCU()->Send();
  receiveRemotes(m);
  bcastRemotes(m,senders);
}

void moveEntities(
    Mesh2* m,
    EntityVector senders[4])
{
  DynamicArray<MeshTag*> tags;
  m->getTags(tags);
  int maxDimension = m->getDimension();
  for (int dimension = 0; dimension <= maxDimension; ++dimension)
  {
    m->getPCU()->Begin();
    sendEntities(m,senders[dimension],tags);
    m->getPCU()->Send();
    EntityVector received;
    receiveEntities(m,tags,received);
    setupRemotes(m,received,senders[dimension]);
  }
}

/* before this call senders are matched to one another
   an no one else, and they each have the correct remote copies
   for their abstract entity.
   Get the full set of copies from each sender, unite between all
   senders, and the result is the new set of matched copies.
   at the end of this routine all senders have correct matching */
static void updateSenderMatching(
    Mesh2* m,
    EntityVector affected[4],
    EntityVector senders[4])
{
  m->getPCU()->Begin();
  int self = m->getPCU()->Self();
  for (int d=0; d < 4; ++d)
  {
    for (size_t i=0; i < senders[d].size(); ++i)
    {
      MeshEntity* e = senders[d][i];
      Matches matches;
      m->getMatches(e,matches);
      if ( ! matches.getSize()) continue;
      Copies copies;
      m->getRemotes(e,copies);
      Parts residence;
      m->getResidence(e,residence);
      /* pack the remote copies to itself, then
         add itself to the copies if it remains */
      m->getPCU()->Pack(self,e);
      packCopies(self,copies, m->getPCU());
      if (residence.count(self))
        copies[self]=e;
      for (size_t j=0; j < matches.getSize(); ++j)
      { /* pack all copies to matched senders */
        m->getPCU()->Pack(matches[j].peer,matches[j].entity);
        packCopies(matches[j].peer,copies, m->getPCU());
      }
    }
/* destroy all matching in the entire affected set.
   It will be restored for senders in the following
   receive and everywhere else by bcastMatching */
    for (size_t i=0; i < affected[d].size(); ++i)
      m->clearMatches(affected[d][i]);
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive())
  {
    MeshEntity* e;
    m->getPCU()->Unpack(e);
    Copies copies;
    unpackCopies(copies, m->getPCU());
    APF_ITERATE(Copies,copies,cit)
      m->addMatch(e,cit->first,cit->second);
  }
}

/* now senders just broadcast the full matching to their remote copies */
static void bcastMatching(
    Mesh2* m,
    EntityVector senders[4])
{
  m->getPCU()->Begin();
  int self = m->getPCU()->Self();
  for (int d=0; d < 4; ++d)
  {
    for (size_t i=0; i < senders[d].size(); ++i)
    {
      MeshEntity* e = senders[d][i];
      Matches matches;
      m->getMatches(e,matches);
      if ( ! matches.getSize()) continue;
      Copies remotes;
      m->getRemotes(e,remotes);
      Parts residence;
      m->getResidence(e,residence);
      if (residence.count(self))
      {
        Match selfMatch;
        selfMatch.peer = self;
        selfMatch.entity = e;
        matches.append(selfMatch);
      }
      APF_ITERATE(Copies,remotes,rit)
      {
        int to = rit->first;
        m->getPCU()->Pack(to,rit->second);
        size_t n = matches.getSize();
        m->getPCU()->Pack(to,n);
        for (size_t j=0; j < n; ++j)
        {
          m->getPCU()->Pack(to,matches[j].peer);
          m->getPCU()->Pack(to,matches[j].entity);
        }
      }
    }
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive())
  {
    MeshEntity* e;
    m->getPCU()->Unpack(e);
    Match match;
    size_t n;
    m->getPCU()->Unpack(n);
    for (size_t i=0; i < n; ++i)
    {
      m->getPCU()->Unpack(match.peer);
      m->getPCU()->Unpack(match.entity);
      if ( ! ((match.peer == self)&&
              (match.entity == e)))
        m->addMatch(e,match.peer,match.entity);
    }
  }
}

/* matched senders now have updated remote copies,
   we unite these sets among all matched senders,
   including the senders themselves if they are
   not being removed, and then broadcast from the
   senders again, this time with matchings */
void updateMatching(
    Mesh2* m,
    EntityVector affected[4],
    EntityVector senders[4])
{
  if ( ! m->hasMatching()) return;
  updateSenderMatching(m,affected,senders);
  bcastMatching(m,senders);
}

void deleteOldEntities(
    Mesh2* m,
    EntityVector affected[4])
{
  int rank = m->getPCU()->Self();
  int maxDimension = m->getDimension();
  for (int d=maxDimension; d >= 0; --d)
    APF_ITERATE(EntityVector,affected[d],it)
    {
      MeshEntity* e = *it;
      Parts residence;
      m->getResidence(e,residence);
      if ( ! residence.count(rank))
        m->destroy(e);
    }
}

/* this is the main migration routine */
static void migrate1(Mesh2* m, Migration* plan)
{
  EntityVector affected[4];
  getAffected(m,plan,affected);
  EntityVector senders[4];
  getSenders(m,affected,senders);
  reduceMatchingToSenders(m,senders);
  updateResidences(m,plan,affected);
  delete plan;
  moveEntities(m,senders);
  updateMatching(m,affected,senders);
  deleteOldEntities(m,affected);
  m->acceptChanges();
}

const size_t maxMigrationLimit = 10*1000*1000;
static size_t migrationLimit = maxMigrationLimit;

void setMigrationLimit(size_t maxElements, pcu::PCU *PCUObj)
{
  if( maxElements >= maxMigrationLimit ) {
    if(!PCUObj->Self())
      lion_eprint(1, "ERROR requested migration limit exceeds"
                      " %lu... exiting\n", maxMigrationLimit);
    abort();
  }
  migrationLimit = maxElements;
}

/* this implements partial migrations
   to limit peak memory use */
static void migrate2(Mesh2* m, Migration* plan)
{
  std::vector<std::pair<MeshEntity*, int> > tmp;
  tmp.resize(plan->count());
  for (size_t i = 0; i < tmp.size(); ++i)
  {
    MeshEntity* e = plan->get(i);
    tmp[i].first = e;
    tmp[i].second = plan->sending(e);
  }
  delete plan;
  size_t sent = 0;
  while (m->getPCU()->Or(sent < tmp.size()))
  {
    plan = new Migration(m);
    size_t send = std::min(tmp.size() - sent, migrationLimit);
    for (size_t i = sent; i < sent + send; ++i)
      plan->send(tmp[i].first, tmp[i].second);
    migrate1(m, plan);
    sent += send;
  }
}

void migrateSilent(Mesh2* m, Migration* plan)
{
  if (m->getPCU()->Or(static_cast<size_t>(plan->count()) > migrationLimit))
    migrate2(m, plan);
  else
    migrate1(m, plan);
}

void migrate(Mesh2* m, Migration* plan)
{
  migrateSilent(m, plan);
  warnAboutEmptyParts(m);
}

}//namespace apf
