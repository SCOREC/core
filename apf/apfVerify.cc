#include "apfMesh.h"
#include "apf.h"
#include <gmi.h>
#include <sstream>
#include <apfGeometry.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include "stdlib.h" // malloc

namespace apf {

static void intersect(
    std::set<int>& a,
    std::set<int> const& b)
{
  for (std::set<int>::iterator it = a.begin();
      it != a.end();)
  {
    if ( ! b.count(*it))
      a.erase(*(it++));
    else
      ++it;
  }
}

static Parts getCandidateParts(Mesh* m, MeshEntity* e)
{
  Downward d;
  int nd = m->getDownward(e, getDimension(m, e) - 1, d);
  Parts c;
  m->getResidence(d[0], c);
  for (int i = 1; i < nd; ++i)
  {
    Parts dr;
    m->getResidence(d[i], dr);
    intersect(c, dr);
  }
  return c;
}

static bool isSubset(Parts const& a, Parts const& b)
{
  APF_ITERATE(Parts, a, it)
    if (!b.count(*it))
      return false;
  return true;
}

static void verifyDown(Mesh* m, MeshEntity* e, int gd, int ed)
{
  Downward d;
  int nd = m->getDownward(e, ed - 1, d);
  for (int i = 0; i < nd; ++i)
  {
    int dgd = m->getModelType(m->toModel(d[i]));
    PCU_ALWAYS_ASSERT(dgd <= gd);
  }
  Parts r;
  m->getResidence(e, r);
  PCU_ALWAYS_ASSERT(isSubset(r, getCandidateParts(m, e)));
}

typedef std::map<ModelEntity*, int> UpwardCounts;
typedef std::map<ModelEntity*, bool> SideManifoldness;

static void getUpwardCounts(gmi_model* gm, int meshDimension, UpwardCounts& uc)
{
  for (int d = 0; d < meshDimension; ++d) {
    gmi_iter* it = gmi_begin(gm, d);
    gmi_ent* ge;
    while ((ge = gmi_next(gm, it))) {
      gmi_set* up = gmi_adjacent(gm, ge, d + 1);
      uc[(ModelEntity*)ge] = up->n;
      gmi_free_set(up);
    }
    gmi_end(gm, it);
  }
}

static void verifyUp(Mesh* m, UpwardCounts& guc,
    MeshEntity* e, bool abort_on_error)
{
  apf::Up up;
  m->getUp(e, up);
  int upwardCount = up.n;

  /* check for duplicates in the upward list */
  for (int i = 0; i < upwardCount; ++i)
    for (int j = i + 1; j < upwardCount; ++j)
      PCU_ALWAYS_ASSERT(up.e[i] != up.e[j]);

  bool adjacentToUpwardGhost=false;
  for (int i=0; i<upwardCount; ++i)
    if (m->isGhost(up.e[i])) adjacentToUpwardGhost=true;

  int meshDimension = m->getDimension();
  int entityDimension = getDimension(m, e);
  int difference = meshDimension - entityDimension;
  PCU_ALWAYS_ASSERT(difference);
  ModelEntity* ge = m->toModel(e);
  int modelDimension = m->getModelType(ge);
  int modelUpwardCount = guc[ge];
  bool isOnNonManifoldFace = (modelDimension == meshDimension - 1) &&
                             (modelUpwardCount > 1);
  bool isOnManifoldBoundary = ( ! isOnNonManifoldFace) &&
                              (modelDimension < meshDimension);
  bool isShared = m->isShared(e);
  bool isMatched = false;
  if (m->hasMatching()) {
    apf::Matches matches;
    m->getMatches(e, matches);
    isMatched = (matches.getSize() > 0);
  }
  bool isExposedByModel = isOnManifoldBoundary;
  bool isExposedByMesh = isShared || isMatched;
  bool isExposed = isExposedByModel || isExposedByMesh;
  bool isOnEqualOrder = (modelDimension == entityDimension);
  int expected;
  if (isExposed)
    expected = difference;
  else
    expected = difference + 1;
  if (!isExposedByMesh && isOnEqualOrder)
    expected = std::max(expected, modelUpwardCount);
  bool okay;
  if (difference == 1)
    okay = (upwardCount == expected);
  else
    okay = (upwardCount >= expected);
  if ( ! okay) 
  {
    std::stringstream ss;
    char const* n = apf::Mesh::typeName[m->getType(e)];
    ss << "apf::Verify: " << n << " with " << upwardCount << " adjacent "
       << dimName[entityDimension + 1] << "s\n";
    ss << "centroid: " << apf::getLinearCentroid(m, e) << '\n';
    ss << "based on the following:\n";
    ss << " - "  << n << " is classified on a model "
       << dimName[modelDimension] << '\n';
    if (modelDimension < meshDimension) {
      ss << "   "  << "which is adjacent to " << modelUpwardCount
         << " model " << dimName[modelDimension + 1] << "s\n";
    }
    if (isOnNonManifoldFace)
      ss << "   " << "and is a non-manifold " << dimName[modelDimension] << '\n';
    if (isOnManifoldBoundary)
      ss << "   " << "and is a manifold boundary " << dimName[modelDimension] << '\n';
    if (isExposedByModel)
      ss << "   making the " << n << " \"exposed\" at the geometric boundary\n";
    if (isShared)
      ss << " - " << n << " is shared in parallel\n";
    if (isMatched)
      ss << " - " << n << " has periodic matches\n";
    if (isExposedByMesh)
      ss << "   making it \"exposed\" at the part boundary\n";
    if (m->isGhost(e))
      ss << " - " << n << " is a ghost entity\n";
    if (adjacentToUpwardGhost)
      ss << " - " << n << " is adjacent to "<<entityDimension+1<<"-dimensional ghost entity\n";
    else
    {
      ss << "we would expect the adjacent " << dimName[entityDimension + 1] << " count to be ";
      if (difference == 1)
        ss << "exactly " << expected << '\n';
      else
        ss << "at least " << expected << '\n';
    }
    std::string s = ss.str(); 
    if (!adjacentToUpwardGhost && abort_on_error)
      fail(s.c_str()); 
  }
  /* this is here for some spiderwebby simmetrix meshes */
  if (upwardCount >= 200) {
    std::stringstream ss;
    ss << "warning: entity of type " << m->getType(e)
      << " at " << getLinearCentroid(m, e) << " has "
      << upwardCount << " upward adjacencies\n";
    // use a stringstream to prevent output from different procs mixing
    std::string s = ss.str();
    lion_eprint(1,"%s",s.c_str());
  }
}

static void verifyResidence(Mesh* m, MeshEntity* e)
{
  Parts p;
  m->getResidence(e, p);
  PCU_ALWAYS_ASSERT(p.count(m->getOwner(e)));
  Copies r;
  m->getRemotes(e, r);
  PCU_ALWAYS_ASSERT(r.size() + 1 == p.size());
  APF_ITERATE(Copies, r, it)
    PCU_ALWAYS_ASSERT(p.count(it->first));
}

static void verifyEntity(Mesh* m, UpwardCounts& guc, MeshEntity* e, bool abort_on_error)
{
  int ed = getDimension(m, e);
  int md = m->getDimension();
  PCU_ALWAYS_ASSERT(md >= ed);
  int gd = m->getModelType(m->toModel(e));
  PCU_ALWAYS_ASSERT(gd >= ed);
  if (ed)
    verifyDown(m, e, gd, ed);
  if (ed < md)
    verifyUp(m, guc, e, abort_on_error);
  verifyResidence(m, e);
}

static Copies getAllCopies(Mesh* m, MeshEntity* e)
{
  Copies c;
  m->getRemotes(e, c);
  c[m->getPCU()->Self()] = e;
  return c;
}

static void verifyAllCopies(Copies& a, pcu::PCU *PCUObj)
{
  APF_ITERATE(Copies, a, it)
  {
    PCU_ALWAYS_ASSERT(it->first >= 0);
    PCU_ALWAYS_ASSERT(it->first < PCUObj->Peers());
  }
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

static void sendAllCopies(Mesh* m, MeshEntity* e)
{
  Copies a;
  Copies r;
  a = getAllCopies(m, e);
  verifyAllCopies(a, m->getPCU());
  m->getRemotes(e, r);
  PCU_ALWAYS_ASSERT(!r.count(m->getPCU()->Self()));
  APF_ITERATE(Copies, r, it)
  {
    m->getPCU()->Pack(it->first, it->second);
    packCopies(it->first, a, m->getPCU());
  }
}

static void receiveAllCopies(Mesh* m)
{
  MeshEntity* e;
  m->getPCU()->Unpack(e);
  Copies a;
  unpackCopies(a, m->getPCU());
  Copies b = getAllCopies(m, e);
  PCU_ALWAYS_ASSERT(a == b);
}

static void verifyRemoteCopies(Mesh* m)
{
  for (int d = 0; d <= m->getDimension(); ++d)
  {
    m->getPCU()->Begin();
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    while ((e = m->iterate(it)))
      if (m->isShared(e) && !m->isGhost(e))
        sendAllCopies(m, e);
    m->end(it);
    m->getPCU()->Send();
    while (m->getPCU()->Receive())
     receiveAllCopies(m);
  }
}

// ghost verification
static void sendGhostCopies(Mesh* m, MeshEntity* e)
{
  Copies g;
  m->getGhosts(e, g);
  PCU_ALWAYS_ASSERT(g.size());
  APF_ITERATE(Copies, g, it)
  {
    m->getPCU()->Pack(it->first, it->second);
    m->getPCU()->Pack(it->first, e);
  }
}

static void receiveGhostCopies(Mesh* m)
{
  int from = m->getPCU()->Sender();
  MeshEntity* e;
  m->getPCU()->Unpack(e);
  MeshEntity* g;
  m->getPCU()->Unpack(g);
  PCU_ALWAYS_ASSERT(m->isGhost(e) || m->isGhosted(e));
  Copies ghosts;
  m->getGhosts(e,ghosts);
  if (m->isGhosted(e))
  {
    PCU_ALWAYS_ASSERT(ghosts.count(from));
  }
  if (m->isGhost(e) && m->getOwner(e)==from)
  {
    PCU_ALWAYS_ASSERT(ghosts.size()==1);
  }
}

static void verifyGhostCopies(Mesh* m)
{
  m->getPCU()->Begin();
  for (int d = 0; d <= m->getDimension(); ++d)
  {
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    while ((e = m->iterate(it)))
      if (m->isGhosted(e) || m->isGhost(e))
        sendGhostCopies(m, e);
    m->end(it);
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive())
    receiveGhostCopies(m);
}

static bool hasMatch(
    Matches& matches,
    int peer,
    MeshEntity* entity)
{
  unsigned i;
  for (i = 0; i < matches.getSize(); ++i)
    if (matches[i].peer == peer &&
        matches[i].entity == entity)
      return true;
  return false;
}

static void sendSelfToMatches(MeshEntity* e, Matches& matches, pcu::PCU *PCUObj)
{
  APF_ITERATE(Matches, matches, it)
  {
    PCU_ALWAYS_ASSERT(!((it->peer == PCUObj->Self())&&(it->entity == e)));
    PCUObj->Pack(it->peer, e);
    PCUObj->Pack(it->peer, it->entity);
  }
}

static void receiveMatches(Mesh* m)
{
  MeshEntity* source;
  m->getPCU()->Unpack(source);
  MeshEntity* e;
  m->getPCU()->Unpack(e);
  Matches matches;
  m->getMatches(e, matches);
  PCU_ALWAYS_ASSERT(hasMatch(matches, m->getPCU()->Sender(), source));
}

static bool hasDuplicates(Matches const& matches) {
  for (size_t i = 0; i < matches.getSize(); ++i)
  for (size_t j = 0; j < matches.getSize(); ++j)
    if (i != j && matches[i].peer == matches[j].peer &&
                  matches[i].entity == matches[j].entity)
      return true;
  return false;
}

static void verifyMatches(Mesh* m)
{
  m->getPCU()->Begin();
  for (int d = 0; d <= m->getDimension(); ++d)
  {
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    while ((e = m->iterate(it))) {
      Matches matches;
      m->getMatches(e, matches);
      PCU_ALWAYS_ASSERT(!hasDuplicates(matches));
      sendSelfToMatches(e, matches, m->getPCU());
    }
    m->end(it);
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive())
    receiveMatches(m);
}

static void sendCoords(Mesh* m, MeshEntity* e)
{
  Vector3 x;
  m->getPoint(e, 0, x);
  Vector3 p(0,0,0);
  m->getParam(e, p);
  Copies r;
  m->getRemotes(e, r);
  APF_ITERATE(Copies, r, it)
  {
    m->getPCU()->Pack(it->first, it->second);
    m->getPCU()->Pack(it->first, x);
    m->getPCU()->Pack(it->first, p);
  }
}

static bool receiveCoords(Mesh* m)
{
  MeshEntity* e;
  Vector3 ox;
  Vector3 op;
  m->getPCU()->Unpack(e);
  m->getPCU()->Unpack(ox);
  m->getPCU()->Unpack(op);
  Vector3 x;
  Vector3 p(0,0,0);
  m->getPoint(e, 0, x);
  m->getParam(e, p);
  return areClose(x, ox, 0.0) &&
         areClose(p, op, 0.0);
}

static long verifyCoords(Mesh* m)
{
  m->getPCU()->Begin();
  MeshIterator* it = m->begin(0);
  MeshEntity* e;
  while ((e = m->iterate(it)))
    if (m->isShared(e))
      sendCoords(m, e);
  m->end(it);
  m->getPCU()->Send();
  long n = 0;
  while (m->getPCU()->Receive())
    if (!receiveCoords(m))
      ++n;
  return m->getPCU()->Add<long>(n);
}

long verifyVolumes(Mesh* m, bool printVolumes)
{
  MeshIterator* it = m->begin(m->getDimension());
  MeshEntity* e;
  long n = 0;
  while ((e = m->iterate(it)))
  {
    if (!isSimplex(m->getType(e)))
      continue;
    double v = measure(m,e);
    if (v < 0) {
      if (printVolumes) {
        std::stringstream ss;
        ss << "warning: element volume " << v
          << " at " << getLinearCentroid(m, e) << '\n';
        std::string s = ss.str();
        lion_oprint(1, "%s", s.c_str());
        fflush(stdout);
      }
      ++n;
    }
  }
  m->end(it);
  return m->getPCU()->Add<long>(n);
}

static void packAlignment(Mesh* m, MeshEntity* e, MeshEntity* r, int to)
{
  m->getPCU()->Pack(to,r);
  int d = getDimension(m, e);
  Downward down;
  int nd = m->getDownward(e, d - 1, down);
  for (int i = 0; i < nd; ++i) {
    Copies remotes;
    m->getRemotes(down[i], remotes);
    MeshEntity* dr = remotes[to];
    m->getPCU()->Pack(to,dr);
  }
}

static void sendAlignment(Mesh* m, MeshEntity* e)
{
  Copies remotes;
  m->getRemotes(e, remotes);
  APF_ITERATE(Copies, remotes, it)
    packAlignment(m, e, it->second, it->first);
}

static void receiveAlignment(Mesh* m)
{
  MeshEntity* e;
  m->getPCU()->Unpack(e);
  int d = getDimension(m, e);
  Downward down;
  int nd = m->getDownward(e, d - 1, down);
  for (int i = 0; i < nd; ++i) {
    m->getPCU()->Unpack(e);
    PCU_ALWAYS_ASSERT(down[i] == e);
  }
}

static void verifyAlignment(Mesh* m)
{
  for (int d = 1; d <= m->getDimension(); ++d) {
    m->getPCU()->Begin();
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    while ((e = m->iterate(it)))
      if (m->isShared(e))
        sendAlignment(m, e);
    m->end(it);
    m->getPCU()->Send();
    while (m->getPCU()->Receive())
      receiveAlignment(m);
  }
}

void packFieldInfo(Field* f, int to, pcu::PCU *PCUObj)
{
  std::string name; 
  name = getName(f);
  packString(name, to, PCUObj);
  int type, size;
  type = getValueType(f);
  PCUObj->Pack(to, type);
  size = countComponents(f);
  PCUObj->Pack(to, size);
}

void unpackFieldInfo(std::string& name, int& type, int& size, pcu::PCU *PCUObj)
{
  name = unpackString(PCUObj);
  PCUObj->Unpack(type);
  PCUObj->Unpack(size);
}

static void verifyFields(Mesh* m) 
{
  int self = m->getPCU()->Self();
  int n = m->countFields();
  std::vector<Field*> fields;
  for (int i=0; i<n; ++i)
    fields.push_back(m->getField(i));
  if (!fields.size()) return;

  m->getPCU()->Begin();
  if (self) {
    m->getPCU()->Pack(self - 1, n);
    for (int i = 0; i < n; ++i)
      packFieldInfo(fields[i], self - 1, m->getPCU());
  }

  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    int n;
    m->getPCU()->Unpack(n);
    PCU_ALWAYS_ASSERT(fields.size() == (size_t)n);
    for (int i = 0; i < n; ++i) {
      std::string name;
      int type;
      int size;
      unpackTagInfo(name, type, size, m->getPCU());
      PCU_ALWAYS_ASSERT(name == getName(fields[i]));
      PCU_ALWAYS_ASSERT(type == getValueType(fields[i]));
      PCU_ALWAYS_ASSERT(size == countComponents(fields[i]));
    }
  }
}

// VERIFY TAGS
static void sendTagData(Mesh* m, MeshEntity* e, DynamicArray<MeshTag*>& tags, Copies& copies, bool is_ghost=false)
{
  void* msg_send;
  MeshEntity** s_ent;
  size_t msg_size;
  int n = tags.getSize(), tag_type, tag_size;

  for (int i = 0; i < n; ++i)
  {
    if (m->getTagName(tags[i])==std::string("ghosted_tag") && is_ghost) continue;
    if (m->getTagName(tags[i])==std::string("ghost_tag")) continue;
    if (!m->hasTag(e, tags[i])) continue;

    APF_ITERATE(Copies, copies, it)
    {
      tag_type = m->getTagType(tags[i]);
      tag_size = m->getTagSize(tags[i]);
      switch (tag_type)
      {
        case Mesh::DOUBLE:  msg_size=sizeof(MeshEntity*)+sizeof(int)+tag_size*sizeof(double); break;
        case Mesh::INT:  msg_size=sizeof(MeshEntity*)+sizeof(int)+tag_size*sizeof(int); break;
        case Mesh::LONG:  msg_size=sizeof(MeshEntity*)+sizeof(int)+tag_size*sizeof(long); break;
        default: msg_size=0;
      }
      PCU_ALWAYS_ASSERT(msg_size);
      msg_send = malloc(msg_size);
      s_ent = (MeshEntity**)msg_send; 
      *s_ent = it->second; 
      int *s_tagid = (int*)((char*)msg_send + sizeof(MeshEntity*));
      s_tagid[0] =  i;
      switch (tag_type)
      {
        case Mesh::DOUBLE: {
                             double *data = (double*)((char*)msg_send+sizeof(MeshEntity*)+sizeof(int));
                             m->getDoubleTag(e, tags[i], data);
                             break;
                           }
        case Mesh::INT:    {
                             int *data = (int*)((char*)msg_send+sizeof(MeshEntity*)+sizeof(int));
                             m->getIntTag(e, tags[i], data);
                             break;
                           }
        case Mesh::LONG:   {
                             long *data = (long*)((char*)msg_send+sizeof(MeshEntity*)+sizeof(int));
                             m->getLongTag(e, tags[i], data);
                             break;
                           }
        default: break;
      }
      m->getPCU()->Write(it->first, (void*)msg_send, msg_size);
      free(msg_send);
    } // apf_iterate
  } // for
}

static void receiveTagData(Mesh* m, DynamicArray<MeshTag*>& tags)
{
  MeshEntity* e;
  int tag_type, tag_size;
  void *msg_recv;
  int pid_from;
  MeshTag* tag;
  size_t msg_size;
  std::set<MeshTag*> mismatch_tags;
  
  while(m->getPCU()->Read(&pid_from, &msg_recv, &msg_size))
  {
    e = *((MeshEntity**)msg_recv); 
    int *id = (int*)((char*)msg_recv+sizeof(MeshEntity*)); 
    tag = tags[id[0]];
    if (mismatch_tags.find(tag)!=mismatch_tags.end()) continue;
    tag_type = m->getTagType(tag);
    tag_size = m->getTagSize(tag);
    if (!m->hasTag(e,tag)) 
    {
      mismatch_tags.insert(tag);
      continue;
    }
    switch (tag_type)
    {
      case Mesh::DOUBLE: {
                           double* r_data = (double*)((char*)msg_recv+sizeof(MeshEntity*)+sizeof(int)); 
                           int num_data = (msg_size-sizeof(MeshEntity*)-sizeof(int))/sizeof(double);
                           PCU_ALWAYS_ASSERT(num_data==tag_size);
                           double *tag_data = new double[tag_size];
                           m->getDoubleTag(e, tag, tag_data);
                           for (int i=0; i<num_data; ++i)
                           {
                             if(tag_data[i]!=r_data[i])
                             {
                                mismatch_tags.insert(tag);
                                break;
                             }
                           }
                           delete [] tag_data;
                           break;
                         }
      case Mesh::INT: {
                           int* r_data = (int*)((char*)msg_recv+sizeof(MeshEntity*)+sizeof(int)); 
                           int num_data = (msg_size-sizeof(MeshEntity*)-sizeof(int))/sizeof(int);
                           PCU_ALWAYS_ASSERT(num_data==tag_size);
                           int *tag_data = new int[tag_size];
                           m->getIntTag(e, tag, tag_data);
                           for (int i=0; i<num_data; ++i)
                           {
                             if (tag_data[i]!=r_data[i])
                             {
                                mismatch_tags.insert(tag);
                                break;
                             }
                           }
                           delete [] tag_data;
                           break;
                         }
      case Mesh::LONG: {
                           long* r_data = (long*)((char*)msg_recv+sizeof(MeshEntity*)+sizeof(int)); 
                           int num_data = (msg_size-sizeof(MeshEntity*)-sizeof(int))/sizeof(long);
                           PCU_ALWAYS_ASSERT(num_data==tag_size);
                           long *tag_data = new long[tag_size];
                           m->getLongTag(e, tag, tag_data);
                           for (int i=0; i<num_data; ++i)
                             if (tag_data[i]!=r_data[i])
                             {
                                mismatch_tags.insert(tag);
                                break;
                             }
                           delete [] tag_data;
                           break;
                         }
        default: break;
    } // switch
  } // while

  int global_size = m->getPCU()->Max<int>((int)mismatch_tags.size());
  if (global_size&&!m->getPCU()->Self())
    for (std::set<MeshTag*>::iterator it=mismatch_tags.begin(); it!=mismatch_tags.end(); ++it)
      lion_oprint(1,"  - tag \"%s\" data mismatch over remote/ghost copies\n", m->getTagName(*it));
}

static void verifyTags(Mesh* m)
{
  DynamicArray<MeshTag*> tags;
  m->getTags(tags);
  int self = m->getPCU()->Self();
  int n = tags.getSize();
  if (!n) return;

  m->getPCU()->Begin();
  if (self) {
    m->getPCU()->Pack(self - 1, n);
    for (int i = 0; i < n; ++i)
      packTagInfo(m, tags[i], self - 1);
  }
  else // master
  {
    if (n)
    {
      lion_oprint(1,"  - verifying tags: ");
      for (int i = 0; i < n; ++i)
      {
        lion_oprint(1,"%s", m->getTagName(tags[i]));
        if (i<n-1) lion_oprint(1,", ");
      }
      lion_oprint(1,"\n");
    }
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    int n;
    m->getPCU()->Unpack(n);
    PCU_ALWAYS_ASSERT(tags.getSize() == (size_t)n);
    for (int i = 0; i < n; ++i) {
      std::string name;
      int type;
      int size;
      unpackTagInfo(name, type, size, m->getPCU());
      PCU_ALWAYS_ASSERT(name == m->getTagName(tags[i]));
      PCU_ALWAYS_ASSERT(type == m->getTagType(tags[i]));
      PCU_ALWAYS_ASSERT(size == m->getTagSize(tags[i]));
    }
  }
  
  // verify tag data

  for (int d = 0; d <= m->getDimension(); ++d) 
  {
    m->getPCU()->Begin();
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    while ((e = m->iterate(it))) 
    {
      if (m->getOwner(e)!=m->getPCU()->Self()) continue;
      if (m->isShared(e)) {
        Copies r;
        m->getRemotes(e, r);
        sendTagData(m, e, tags, r);
      }
      if (m->isGhosted(e)) {
        Copies g;
        m->getGhosts(e, g);
        sendTagData(m, e, tags, g, true);
      }
    } // while
    m->end(it);
    m->getPCU()->Send();
    receiveTagData(m, tags);
  } // for
}

void verify(Mesh* m, bool abort_on_error)
{
  double t0 = pcu::Time();
  verifyTags(m);
  verifyFields(m);

  UpwardCounts guc;
  getUpwardCounts(m->getModel(), m->getDimension(), guc);
  /* got to 3 on purpose, so we can verify if
     m->getDimension is lying */
  for (int d = 0; d <= 3; ++d)
  {
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    size_t n = 0;
    while ((e = m->iterate(it)))
    {
      verifyEntity(m, guc, e, abort_on_error);
      ++n;
    }
    m->end(it);
    PCU_ALWAYS_ASSERT(n == m->count(d));
    if (d > m->getDimension())
      PCU_ALWAYS_ASSERT(!n);
  }
  guc.clear();
  verifyRemoteCopies(m);
  verifyGhostCopies(m);
  verifyAlignment(m);
  verifyMatches(m);
  long n = verifyCoords(m);
  if (n && (!m->getPCU()->Self()))
    lion_eprint(1,"apf::verify fail: %ld coordinate mismatches\n", n);
  n = verifyVolumes(m);
  if (n && (!m->getPCU()->Self()))
    lion_eprint(1,"apf::verify warning: %ld negative simplex elements\n", n);
  double t1 = pcu::Time();
  if (!m->getPCU()->Self())
    lion_oprint(1,"mesh verified in %f seconds\n", t1 - t0);
}

}
