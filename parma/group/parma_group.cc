#include <PCU.h>
#include <parma.h>

using apf::Remap;

static void retreat(apf::Mesh2* m, Remap& remap)
{
  int to = remap(PCU_Comm_Self());
  apf::Migration* plan = new apf::Migration(m);
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* e;
  while ((e = m->iterate(it)))
    plan->send(e, to);
  m->end(it);
  apf::migrateSilent(m, plan);
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
  if (m)
    apf::remapPartition(m, inMap);
  code.run(group);
  PCU_Switch_Comm(oldComm);
  MPI_Comm_free(&groupComm);
  if (m)
    apf::remapPartition(m, outMap);
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
  apf::Divide inMap(factor);
  apf::Modulo groupMap(factor);
  apf::Round retreatMap(factor);
  apf::Multiply outMap(factor);
  retreatToGroup(m, factor, inMap, retreatMap, groupMap, outMap, toRun);
}

void Parma_SplitPartition(apf::Mesh2* m, int factor, Parma_GroupCode& toRun)
{
  apf::Modulo inMap(factor);
  apf::Divide groupMap(factor);
  apf::Unmodulo outMap(PCU_Comm_Self(), factor);
  runInGroups(m, inMap, groupMap, outMap, toRun);
}

