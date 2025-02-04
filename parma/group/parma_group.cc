#include <parma.h>
#include <pcu_util.h>
#include <memory>

using apf::Remap;

static void retreat(apf::Mesh2* m, Remap& remap)
{
  int to = remap(m->getPCU()->Self());
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
    pcu::PCU *PCUObj,
    Remap& inMap,
    Remap& groupMap,
    Remap& outMap,
    GroupCode& code
    )
{
  int self;
  pcu::PCU* expandedPCU;
  if(PCUObj == nullptr){
    self = m->getPCU()->Self();
    expandedPCU = m->getPCU();
  } else {
    self = PCUObj->Self();
    expandedPCU = PCUObj;
  }
  
  int groupRank = inMap(self);
  int group = groupMap(self);
  
  MPI_Comm groupComm;
  PCU_Comm_Split(expandedPCU->GetMPIComm(), group, groupRank, &groupComm);
  auto groupedPCU = std::unique_ptr<pcu::PCU>(new pcu::PCU(groupComm));
  if (m){
    m->switchPCU(groupedPCU.get());
    apf::remapPartition(m, inMap);
  }
  code.PCUObj = std::move(groupedPCU);
  code.run(group);
  PCU_Comm_Free_One(&groupComm);
  MPI_Comm_free(&groupComm);
  if (m){
    m->switchPCU(expandedPCU);
    apf::remapPartition(m, outMap);
  }
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
    GroupCode& code,
    pcu::PCU *PCUObj)
{
  retreat(m, retreatMap);
  RetreatCode retreatCode(code, m, factor);
  runInGroups(m, PCUObj, inMap, groupMap, outMap, retreatCode);
  m->migrate(retreatCode.plan);
}

void Parma_ShrinkPartition(apf::Mesh2* m, int factor, Parma_GroupCode& toRun, pcu::PCU *PCUObj)
{
  if(m == nullptr){
    PCU_ALWAYS_ASSERT(PCUObj != nullptr);
  }
  if(m != nullptr && PCUObj != nullptr){
    PCU_ALWAYS_ASSERT(m->getPCU() == PCUObj);
  }
  apf::Divide inMap(factor);
  apf::Modulo groupMap(factor);
  apf::Round retreatMap(factor);
  apf::Multiply outMap(factor);
  retreatToGroup(m, factor, inMap, retreatMap, groupMap, outMap, toRun, PCUObj);
}

void Parma_SplitPartition(apf::Mesh2* m, int factor, Parma_GroupCode& toRun, pcu::PCU *PCUObj)
{
  if(m == nullptr){
    PCU_ALWAYS_ASSERT(PCUObj != nullptr);
  }
  if(m != nullptr && PCUObj != nullptr){
    PCU_ALWAYS_ASSERT(m->getPCU() == PCUObj);
  }
  apf::Modulo inMap(factor);
  apf::Divide groupMap(factor);
  if(PCUObj == nullptr){
    apf::Unmodulo outMap(m->getPCU()->Self(), factor);
    runInGroups(m, PCUObj, inMap, groupMap, outMap, toRun);
  } else {
    apf::Unmodulo outMap(PCUObj->Self(), factor);
    runInGroups(m, PCUObj, inMap, groupMap, outMap, toRun);
  }
}

