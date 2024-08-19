#include "pcu_util.h"
#include "apfConvert.h"
#include "apfMesh2.h"
#include "apf.h"
#include "apfNumbering.h"
#include <map>
#include <lionPrint.h>

namespace apf {

static void constructVerts(
    Mesh2* m, const Gid* conn, int nelem, int etype,
    GlobalToVert& result)
{
  ModelEntity* interior = m->findModelEntity(m->getDimension(), 0);
  int end = nelem * apf::Mesh::adjacentCount[etype][0];
  int self2 = m->getPCU()->Self();
  for (int i = 0; i < end; ++i)
    if ( ! result.count(conn[i])) {
      result[conn[i]] = m->createVert_(interior);
//      if(conn[i] < 0 || conn[i] > 4305368187 ) { // for whatever reason max is not stored but is found and checked later
      if(conn[i] < 0 ) { 
        lion_eprint(1, "constructVerts building globalToVert: self=%d,gid=%ld,i=%d,nelem=%ld  \n",self2,conn[i],i,nelem);
      }
    }
}

static NewElements constructElements(
    Mesh2* m, const Gid* conn, int nelem, int etype,
    GlobalToVert& globalToVert)
{
  ModelEntity* interior = m->findModelEntity(m->getDimension(), 0);
  int nev = apf::Mesh::adjacentCount[etype][0];
  NewElements newElements;
  for (int i = 0; i < nelem; ++i) {
    Downward verts;
    int offset = i * nev;
    for (int j = 0; j < nev; ++j)
      verts[j] = globalToVert[conn[j + offset]];
    newElements.push_back(buildElement(m, interior, etype, verts));
  }
  return newElements;
}

static Gid getMax(const GlobalToVert& globalToVert, Mesh2* m)
{
  Gid max = -1;
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it)
    max = std::max(max, it->first);
  return m->getPCU()->Max<long>(max); // this is type-dependent
}


/* algorithm courtesy of Sebastian Rettenberger:
   use brokers/routers for the vertex global ids.
   Although we have used this trick before (see mpas/apfMPAS.cc),
   I didn't think to use it here, so credit is given. */
static void constructResidence(Mesh2* m, GlobalToVert& globalToVert)
{
  Gid ifirst=0;
  int self2 = m->getPCU()->Self();
  APF_ITERATE(GlobalToVert, globalToVert, it) {
    Gid gid = it->first;  
    if(ifirst==0 || ifirst==13437400 ) {
        lion_eprint(1, "constructResidence: self=%d,gid=%ld,ifirst=%ld  \n",self2,gid,ifirst);
    }
    ifirst++;
  }
  Gid max = getMax(globalToVert, m);  // seems like we read this and know it already on every rank so why compute with global comm?
  ifirst=0;
  APF_ITERATE(GlobalToVert, globalToVert, it) {
    Gid gid = it->first;  
    if(gid < 0 || gid > max ){ 
        lion_eprint(1, "constructResidence cgTV2: self=%d,gid=%ld,ifirst=%ld  \n",self2,gid,ifirst);
    }
    if(ifirst==0 || ifirst==13437400 ) {
        lion_eprint(1, "constructResidence: self=%d,gid=%ld,ifirst=%ld,max=%ld  \n",self2,gid,ifirst,max);
    }
    ifirst++;
  }
  Gid total = max + 1;
  int peers = m->getPCU()->Peers();
  Gid quotientL = total / peers; 
  int quotient = quotientL;
  Gid remainderL = total % peers; 
  int remainder = remainderL;
  int mySize = quotient;
  int self = m->getPCU()->Self();
  if (self == (peers - 1))
    mySize += remainder;
  if (self == (peers - 1)) lion_eprint(1, "CR1 mysize=%d \n",mySize);
  typedef std::vector< std::vector<int> > TmpParts;
  TmpParts tmpParts(mySize);
  /* if we have a vertex, send its global id to the
     broker for that global id */
  m->getPCU()->Begin();
  APF_ITERATE(GlobalToVert, globalToVert, it) {
    Gid gid = it->first;  
    if(gid < 0 || gid > max ){ 
        lion_eprint(1, "constructResidence cgTV3: self=%d,gid=%ld \n",self2,gid);
    } 
    Gid tmpL=gid / quotient;
    int tmpI=tmpL;    int to = std::min(peers - 1,tmpI); 
    m->getPCU()->Pack(to, gid);
  }
  m->getPCU()->Send();
  Gid myOffset = (long)self * quotient;
  if (self == (peers - 1)) lion_eprint(1, "CR5: self=%d,myOffset=%ld,quotient=%d \n",self,myOffset,quotient);
  /* brokers store all the part ids that sent messages
     for each global id */
  while (m->getPCU()->Receive()) {
    Gid gid;
    m->getPCU()->Unpack(gid);
    int from = m->getPCU()->Sender();
    Gid tmpL=gid - myOffset; // forcing 64 bit difference until we know it is safe
    int tmpI=tmpL;
    tmpParts.at(tmpI).push_back(from);
  }
  /* for each global id, send all associated part ids
     to all associated parts */
  m->getPCU()->Begin();
  for (int i = 0; i < mySize; ++i) {
    std::vector<int>& parts = tmpParts[i];
    for (size_t j = 0; j < parts.size(); ++j) {
      int to = parts[j];
      Gid gid = i + myOffset;
      int nparts = parts.size();
      m->getPCU()->Pack(to, gid);
      m->getPCU()->Pack(to, nparts);
      for (size_t k = 0; k < parts.size(); ++k)
        m->getPCU()->Pack(to, parts[k]);
    }
  }
  m->getPCU()->Send();
  /* receiving a global id and associated parts,
     lookup the vertex and classify it on the partition
     model entity for that set of parts */
  while (m->getPCU()->Receive()) {
    Gid gid;
    m->getPCU()->Unpack(gid);
    int nparts;
    m->getPCU()->Unpack(nparts);
    Parts residence;
    for (int i = 0; i < nparts; ++i) {
      int part;
      m->getPCU()->Unpack(part);
      residence.insert(part);
    }
    MeshEntity* vert = globalToVert[gid];
    m->setResidence(vert, residence);
  }
}

/* given correct residence from the above algorithm,
   negotiate remote copies by exchanging (gid,pointer)
   pairs with parts in the residence of the vertex */
static void constructRemotes(Mesh2* m, GlobalToVert& globalToVert)
{
  int self = m->getPCU()->Self();
  m->getPCU()->Begin();
  APF_ITERATE(GlobalToVert, globalToVert, it) {
    Gid gid = it->first;
    if(gid < 0 )      
        lion_eprint(1, "constructRemotes cgTV4: self=%d,gid=%ld \n",self,gid);
    MeshEntity* vert = it->second;
    Parts residence;
    m->getResidence(vert, residence);
    APF_ITERATE(Parts, residence, rit)
      if (*rit != self) {
        m->getPCU()->Pack(*rit, gid);
        m->getPCU()->Pack(*rit, vert);
      }
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    Gid gid;
    m->getPCU()->Unpack(gid);
    MeshEntity* remote;
    m->getPCU()->Unpack(remote);
    int from = m->getPCU()->Sender();
    MeshEntity* vert = globalToVert[gid];
    m->addRemote(vert, from, remote);
  }
  // who is not stuck?
  lion_eprint(1, "%d done inside remotes \n",m->getPCU()->Self());
}

NewElements assemble(Mesh2* m, const Gid* conn, int nelem, int etype,
     GlobalToVert& globalToVert)
{
  constructVerts(m, conn, nelem, etype, globalToVert);
  return constructElements(m, conn, nelem, etype, globalToVert);
}

void finalise(Mesh2* m, GlobalToVert& globalToVert)
{
  constructResidence(m, globalToVert);
  lion_eprint(1, "%d after residence \n",m->getPCU()->Self());
  constructRemotes(m, globalToVert);
  stitchMesh(m);
  m->acceptChanges();
}
 
NewElements construct(Mesh2* m, const Gid* conn, int nelem, int etype,
    GlobalToVert& globalToVert)
{
  const auto& newElements = assemble(m, conn, nelem, etype, globalToVert);
  finalise(m, globalToVert);
  return newElements;
}

void setCoords(Mesh2* m, const double* coords, int nverts,
    GlobalToVert& globalToVert)
{
  Gid max = getMax(globalToVert, m);
  Gid total = max + 1;
  int peers = m->getPCU()->Peers();
  int quotient = total / peers;
  int remainder = total % peers;
  int mySize = quotient;
  int self = m->getPCU()->Self();
  if (self == (peers - 1))
    mySize += remainder;
  Gid myOffset = (long)self * quotient;

  /* Force each peer to have exactly mySize verts.
     This means we might need to send and recv some coords */
  double* c = new double[mySize*3];

  Gid start = m->getPCU()->Exscan(nverts);

  m->getPCU()->Begin();  // the forced 64 bit math below may not be necessary 
  Gid tmpL=start / quotient; 
  int tmpInt=tmpL;
  int to = std::min(peers - 1, tmpInt);
  tmpL=(to+1)*(long)quotient-start; 
  tmpInt=tmpL;
  int n = std::min(tmpInt, nverts);

  while (nverts > 0) {
    m->getPCU()->Pack(to, start);
    m->getPCU()->Pack(to, n);
    m->getPCU()->Pack(to, coords, n*3*sizeof(double));

    nverts -= n;
    start += n;
    coords += n*3;
    to = std::min(peers - 1, to + 1);
    n = std::min(quotient, nverts);
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    m->getPCU()->Unpack(start);
    m->getPCU()->Unpack(n); // |||||| more in-place 64 bit math
    m->getPCU()->Unpack(&c[(start - myOffset) * 3], n*3*sizeof(double));
  }

  /* Tell all the owners of the coords what we need */
  typedef std::vector< std::vector<int> > TmpParts;
  TmpParts tmpParts(mySize);
  m->getPCU()->Begin();
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it) {
    Gid gid = it->first;
    Gid tmpL=gid / quotient;
    int tmpInt=tmpL;
    int to = std::min(peers - 1, tmpInt);
    m->getPCU()->Pack(to, gid);
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    Gid gid;
    m->getPCU()->Unpack(gid);
    Gid from = m->getPCU()->Sender();
    Gid tmpL=gid - myOffset;
    int tmpInt=tmpL;
    tmpParts.at(tmpInt).push_back(from);
  }

  /* Send the coords to everybody who want them */
  m->getPCU()->Begin();
  for (int i = 0; i < mySize; ++i) {
    std::vector<int>& parts = tmpParts[i];
    for (size_t j = 0; j < parts.size(); ++j) {
      int to = parts[j];
      Gid gid = i + myOffset;
      m->getPCU()->Pack(to, gid);
      m->getPCU()->Pack(to, &c[i*3], 3*sizeof(double));
    }
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    Gid gid;
    m->getPCU()->Unpack(gid);
    double v[3];
    m->getPCU()->Unpack(v, sizeof(v));
    Vector3 vv(v);
    m->setPoint(globalToVert[gid], 0, vv);
  }

  delete [] c;
}

void setMatches(Mesh2* m, const Gid* matches, int nverts,
    GlobalToVert& globalToVert)
{
  Gid max = getMax(globalToVert, m);
  Gid total = max + 1;
  int peers = m->getPCU()->Peers();
  int quotient = total / peers;
  int remainder = total % peers;
  int mySize = quotient;
  int self = m->getPCU()->Self();
  if (self == (peers - 1))
    mySize += remainder;
  Gid myOffset = (long)self * quotient;

  /* Force each peer to have exactly mySize verts.
     This means we might need to send and recv some matches */
  Gid* c = new Gid[mySize];
  Gid start = m->getPCU()->Exscan(nverts);

  m->getPCU()->Begin();

  Gid tmpL=start / quotient; 
  int tmpInt=tmpL;
  int to = std::min(peers - 1, tmpInt);
  tmpL=(to+1)*(long)quotient-start; 
  tmpInt=tmpL;
  int n = std::min(tmpInt, nverts);
  while (nverts > 0) {
    m->getPCU()->Pack(to, start);
    m->getPCU()->Pack(to, n);
    m->getPCU()->Pack(to, matches, n*sizeof(Gid));

    nverts -= n;
    start += n;
    matches += n;
    to = std::min(peers - 1, to + 1);
    n = std::min(quotient, nverts);
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    m->getPCU()->Unpack(start);
    m->getPCU()->Unpack(n); /// in-place 64
    m->getPCU()->Unpack(&c[(start - myOffset)], n*sizeof(Gid));
  }


  /* Tell all the owners of the matches what we need */
  typedef std::vector< std::vector<int> > TmpParts;
  TmpParts tmpParts(mySize);
  m->getPCU()->Begin();
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it) {
    Gid gid = it->first;  
    int tmpI=gid / quotient;
    int to = std::min(peers - 1,tmpI); 
    m->getPCU()->Pack(to, gid);
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    Gid gid;
    m->getPCU()->Unpack(gid);
    int from = m->getPCU()->Sender();
    tmpParts.at(gid - myOffset).push_back(from);
  }
 
  MeshTag* matchGidTag = m->createLongTag("matchGids", 1);
  /* Send the matches to everybody who wants them */
  m->getPCU()->Begin();
  for (int i = 0; i < mySize; ++i) {
    std::vector<int>& parts = tmpParts[i];
    for (size_t j = 0; j < parts.size(); ++j) {
      int to = parts[j];
      Gid gid = i + myOffset;
      Gid matchGid = c[i];
      m->getPCU()->Pack(to, gid);
      m->getPCU()->Pack(to, matchGid);
    }
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    Gid gid;
    m->getPCU()->Unpack(gid);
    Gid match;
    m->getPCU()->Unpack(match);
    PCU_ALWAYS_ASSERT(gid != match);
    PCU_ALWAYS_ASSERT(globalToVert.count(gid));
    m->setLongTag(globalToVert[gid], matchGidTag, &match);
  }

  /* Use the 1D partitioning of global ids to distribute the 
   * entity pointers and their owning ranks. Process 0 will hold
   * the entity pointers and owners for mesh vertex gid [0..quotient),
   * process 1 holds gids [quotient..2*quotient), ...
   */
  m->getPCU()->Begin();
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it) {
    MeshEntity* e = it->second;  // KEJ does not follow this
    Gid gid = it->first;
    int tmpI=gid / quotient;
    int to = std::min(peers - 1,tmpI); 
    m->getPCU()->Pack(to, gid);
    m->getPCU()->Pack(to, e);
  }
  m->getPCU()->Send();
  typedef std::pair< int, apf::MeshEntity* > EntOwnerPtrs;
  typedef std::map< Gid, std::vector< EntOwnerPtrs > > GidPtrs;
  GidPtrs gidPtrs;
  while (m->getPCU()->Receive()) {
    Gid gid;
    m->getPCU()->Unpack(gid);
    MeshEntity* vert;
    m->getPCU()->Unpack(vert);
    int owner = m->getPCU()->Sender();
    gidPtrs[gid-myOffset].push_back(EntOwnerPtrs(owner,vert));
  }

  /* Tell the brokers of the matches we need */
  typedef std::pair<Gid,Gid> MatchingPair;
  typedef std::map< MatchingPair, std::vector<Gid> > MatchMap;
  MatchMap matchParts;
  m->getPCU()->Begin();
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it) { //loop over local verts
    Gid gid = it->first;
    Gid matchGid;
    m->getLongTag(it->second, matchGidTag, &matchGid); //get the matched ent gid
    if( matchGid != -1 ) {  // marker for an unmatched vertex
      int tmpI=matchGid / quotient;
      int to = std::min(peers - 1,tmpI);  //broker
      m->getPCU()->Pack(to, gid); // send the local vert gid
      m->getPCU()->Pack(to, matchGid); // and the match gid needed
    }
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    Gid gid;
    m->getPCU()->Unpack(gid); // request from entity gid
    Gid matchGid;
    m->getPCU()->Unpack(matchGid); // requesting matched entity gid 
    MatchingPair mp(gid,matchGid);
    int from = m->getPCU()->Sender();
    matchParts[mp].push_back(from); // store a list of the proceses that need the pair (entity gid, match gid)
  }

  /* Send the match pointer and owner process to everybody
   * who wants them */
  m->getPCU()->Begin();
  APF_ITERATE(MatchMap,matchParts,it) {
    MatchingPair mp = it->first;
    Gid gid = mp.first;
    Gid matchGid = mp.second;
    std::vector<Gid> parts = it->second;
    for(size_t i=0; i<parts.size(); i++) {
      const int to = parts[i];
      m->getPCU()->Pack(to, gid);
      m->getPCU()->Pack(to, matchGid);
      size_t numMatches = gidPtrs[matchGid-myOffset].size();
      m->getPCU()->Pack(to, numMatches);
      for( size_t i=0; i<gidPtrs[matchGid-myOffset].size(); i++) {
        EntOwnerPtrs eop = gidPtrs[matchGid-myOffset][i];
        int owner = eop.first;
        m->getPCU()->Pack(to, owner);
        apf::MeshEntity* ent = eop.second;
        m->getPCU()->Pack(to, ent);
      }
    }
  }
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    Gid gid;
    m->getPCU()->Unpack(gid);
    Gid matchGid;
    m->getPCU()->Unpack(matchGid);
    size_t numMatches;
    m->getPCU()->Unpack(numMatches);
    for(size_t i=0; i<numMatches; i++) {
      int owner;
      m->getPCU()->Unpack(owner);
      MeshEntity* match;
      m->getPCU()->Unpack(match);
      PCU_ALWAYS_ASSERT(globalToVert.count(gid));
      MeshEntity* partner = globalToVert[gid];
      PCU_ALWAYS_ASSERT(! (match == partner && owner == self) );
      m->addMatch(partner, owner, match);
    }
  }

  APF_CONST_ITERATE(GlobalToVert, globalToVert, it) { //loop over local verts
    apf::MeshEntity* left = it->second;
    Gid matchGid;
    m->getLongTag(left, matchGidTag, &matchGid); //get the matched ent gid
    if( matchGid != -1 ) {  // a matched vtx
      apf::Copies copies;
      m->getRemotes(left,copies);
      APF_ITERATE(apf::Copies, copies, cp) {
        int rightPart = cp->first;
        apf::MeshEntity* right = cp->second;
        m->addMatch(left, rightPart, right);
      }
    }
  }

  delete [] c;

  apf::removeTagFromDimension(m, matchGidTag, 0);
  m->destroyTag(matchGidTag);
}

void destruct(Mesh2* m, Gid*& conn, int& nelem, int &etype, int cellDim)
{
  if(cellDim == -1) cellDim = m->getDimension();
  //int dim = m->getDimension();
  nelem = m->count(cellDim);
  conn = 0;
  GlobalNumbering* global = makeGlobal(numberOwnedNodes(m, "apf_destruct"));
  synchronize(global);
  MeshIterator* it = m->begin(cellDim);
  MeshEntity* e;
  int i = 0;
  while ((e = m->iterate(it))) {
    etype = m->getType(e);
    Downward verts;
    int nverts = m->getDownward(e, 0, verts);
    if (!conn)
      conn = new Gid[nelem * nverts];
    for (int j = 0; j < nverts; ++j)
      conn[i++] = getNumber(global, Node(verts[j], 0));
  }
  m->end(it);
  destroyGlobalNumbering(global);
}

void extractCoords(Mesh2* m, double*& coords, int& nverts)
{
  nverts = countOwned(m, 0);
  coords = new double[nverts * 3];

  MeshIterator* it = m->begin(0);
  int i = 0;
  while (MeshEntity* v = m->iterate(it)) {
    if (m->isOwned(v)) {
      Vector3 p;
      m->getPoint(v, 0, p);
      p.toArray(&coords[i*3]);
      i++;
    }
  }
  m->end(it);
}

}
