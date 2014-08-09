#include <parma.h>
#include <PCU.h>

struct Remap
{
  virtual int operator()(int n) = 0;
};

struct Divide : public Remap
{
  Divide(int n):by(n) {}
  int by;
  int operator()(int n) {return n / by;}
};

struct Multiply : public Remap
{
  Multiply(int n):by(n) {}
  int by;
  int operator()(int n) {return n * by;}
};

struct Modulo : public Remap
{
  Modulo(int n):by(n) {}
  int by;
  int operator()(int n) {return n % by;}
};

struct Round : public Remap
{
  Round(int n):factor(n) {}
  int factor;
  int operator()(int n) {return (n / factor) * factor;}
};

static void remapResidence(apf::Mesh2* m, Remap& remap)
{
  for (int d = 0; d <= m->getDimension(); ++d) {
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      apf::Parts residence;
      m->getResidence(e, residence);
      apf::Parts newResidence;
      APF_ITERATE(apf::Parts, residence, rit)
        newResidence.insert( remap(*rit) );
      m->setResidence(e, newResidence);
    }
    m->end(it);
  }
}

static void remapRemotes(apf::Mesh2* m, Remap& remap)
{
  for (int d = 0; d < m->getDimension(); ++d) {
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      if ( ! m->isShared(e))
        continue;
      apf::Copies remotes;
      m->getRemotes(e, remotes);
      apf::Copies newRemotes;
      APF_ITERATE(apf::Copies, remotes, rit)
        newRemotes[ remap(rit->first) ] = rit->second;
      m->setRemotes(e, newRemotes);
    }
    m->end(it);
  }
}

static void remapMatches(apf::Mesh2* m, Remap& remap)
{
  if (!m->hasMatching())
    return;
  for (int d = 0; d < m->getDimension(); ++d) {
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      apf::Matches matches;
      m->getMatches(e, matches);
      if (!matches.getSize())
        continue;
      m->clearMatches(e);
      for (size_t i = 0; i < matches.getSize(); ++i)
        m->addMatch(e, remap( matches[i].peer ), matches[i].entity);
    }
    m->end(it);
  }
}

static void remapPartition(apf::Mesh2* m, Remap& remap)
{
  remapResidence(m, remap);
  remapRemotes(m, remap);
  remapMatches(m, remap);
  m->acceptChanges();
}

static void retreat(apf::Mesh* m, Remap& remap)
{
  int to = remap(PCU_Comm_Self());
  apf::Migration* plan = new apf::Migration(m);
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* e;
  while ((e = m->iterate(it)))
    plan->send(e, to);
  m->end(it);
  m->migrate(plan);
}

static apf::Migration* planExpansion(apf::Mesh* m, int factor)
{
  apf::Splitter* s = Parma_MakeRibSplitter(m);
  apf::Migration* plan = s->split(NULL, 1.10, factor);
  delete s;
  return plan;
}

typedef Parma_GroupCode GroupCode;

static void runInGroups(
    apf::Mesh2* m,
    Remap& inMap,
    Remap& groupMap,
    Remap& outMap,
    GroupCode& code)
{
  int self = PCU_Comm_Self();
  int groupRank = inMap(self);
  int group = groupMap(self);
  MPI_Comm oldComm = PCU_Get_Comm();
  MPI_Comm groupComm;
  MPI_Comm_split(oldComm, group, groupRank, &groupComm);
  PCU_Switch_Comm(groupComm);
  remapPartition(m, inMap);
  code.run(group);
  PCU_Switch_Comm(oldComm);
  MPI_Comm_free(&groupComm);
  remapPartition(m, outMap);
}

struct RetreatCode : public GroupCode
{
  GroupCode* next;
  apf::Migration* plan;
  apf::Mesh* mesh;
  int factor;
  RetreatCode(GroupCode& c, apf::Mesh* m, int f)
  {
    next = &c;
    plan = 0;
    mesh = m;
    factor = f;
  }
  void run(int group)
  {
    if (!group) {
      next->run(group);
      plan = planExpansion(mesh, factor);
    } else {
      plan = new apf::Migration(mesh);
    }
  }
};

static void retreatToGroup(
    apf::Mesh2* m,
    int factor,
    Remap& inMap,
    Remap& retreatMap,
    Remap& groupMap,
    Remap& outMap,
    GroupCode& code)
{
  retreat(m, retreatMap);
  RetreatCode retreatCode(code, m, factor);
  runInGroups(m, inMap, groupMap, outMap, retreatCode);
  m->migrate(retreatCode.plan);
}

void Parma_ShrinkPartition(apf::Mesh2* m, int factor, Parma_GroupCode& toRun)
{
  Divide inMap(factor);
  Modulo groupMap(factor);
  Round retreatMap(factor);
  Multiply outMap(factor);
  retreatToGroup(m, factor, inMap, retreatMap, groupMap, outMap, toRun);
}

