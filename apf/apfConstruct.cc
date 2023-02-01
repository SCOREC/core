#include <PCU.h>
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
  int self2 = PCU_Comm_Self();
  for (int i = 0; i < end; ++i)
    if ( ! result.count(conn[i])) {
      result[conn[i]] = m->createVert_(interior);
//      if(conn[i] < 0 || conn[i] > 4305368187 ) { // for whatever reason max is not stored but is found and checked later
      if(conn[i] < 0 ) { 
        lion_eprint(1, "constructVerts building globalToVert: self=%d,gid=%ld,i=%d,nelem=%ld  \n",self2,conn[i],i,nelem);
      }
    }
}

static void constructElements(
    Mesh2* m, const Gid* conn, int nelem, int etype,
    GlobalToVert& globalToVert)
{
  ModelEntity* interior = m->findModelEntity(m->getDimension(), 0);
  int nev = apf::Mesh::adjacentCount[etype][0];
  for (int i = 0; i < nelem; ++i) {
    Downward verts;
    int offset = i * nev;
    for (int j = 0; j < nev; ++j)
      verts[j] = globalToVert[conn[j + offset]];
    buildElement(m, interior, etype, verts);
  }
}

static Gid getMax(const GlobalToVert& globalToVert)
{
  Gid max = -1;
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it)
    max = std::max(max, it->first);
  return PCU_Max_Long(max); // this is type-dependent
}


/* algorithm courtesy of Sebastian Rettenberger:
   use brokers/routers for the vertex global ids.
   Although we have used this trick before (see mpas/apfMPAS.cc),
   I didn't think to use it here, so credit is given. */
static void constructResidence(Mesh2* m, GlobalToVert& globalToVert)
{
  Gid ifirst=0;
  int self2 = PCU_Comm_Self();
  APF_ITERATE(GlobalToVert, globalToVert, it) {
    Gid gid = it->first;  
    if(ifirst==0 || ifirst==13437400 ) {
        lion_eprint(1, "constructResidence: self=%d,gid=%ld,ifirst=%ld  \n",self2,gid,ifirst);
    }
    ifirst++;
  }
  PCU_Barrier();
  Gid max = getMax(globalToVert);  // seems like we read this and know it already on every rank so why compute with global comm?
  PCU_Barrier();
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
  int peers = PCU_Comm_Peers();
  Gid quotientL = total / peers; 
  int quotient = quotientL;
  Gid remainderL = total % peers; 
  int remainder = remainderL;
  int mySize = quotient;
  int self = PCU_Comm_Self();
  if (self == (peers - 1))
    mySize += remainder;
  if (self == (peers - 1)) lion_eprint(1, "CR1 mysize=%d \n",mySize);
  typedef std::vector< std::vector<int> > TmpParts;
  TmpParts tmpParts(mySize);
  /* if we have a vertex, send its global id to the
     broker for that global id */
  PCU_Comm_Begin();
  APF_ITERATE(GlobalToVert, globalToVert, it) {
    Gid gid = it->first;  
    if(gid < 0 || gid > max ){ 
        lion_eprint(1, "constructResidence cgTV3: self=%d,gid=%ld \n",self2,gid);
    } 
    Gid tmpL=gid / quotient;
    int tmpI=tmpL;    int to = std::min(peers - 1,tmpI); 
    PCU_COMM_PACK(to, gid);
  }
  PCU_Comm_Send();
  Gid myOffset = (long)self * quotient;
  if (self == (peers - 1)) lion_eprint(1, "CR5: self=%d,myOffset=%ld,quotient=%d \n",self,myOffset,quotient);
  /* brokers store all the part ids that sent messages
     for each global id */
  while (PCU_Comm_Receive()) {
    Gid gid;
    PCU_COMM_UNPACK(gid);
    int from = PCU_Comm_Sender();
    Gid tmpL=gid - myOffset; // forcing 64 bit difference until we know it is safe
    int tmpI=tmpL;
    tmpParts.at(tmpI).push_back(from);
  }
  /* for each global id, send all associated part ids
     to all associated parts */
  PCU_Comm_Begin();
  for (int i = 0; i < mySize; ++i) {
    std::vector<int>& parts = tmpParts[i];
    for (size_t j = 0; j < parts.size(); ++j) {
      int to = parts[j];
      Gid gid = i + myOffset;
      int nparts = parts.size();
      PCU_COMM_PACK(to, gid);
      PCU_COMM_PACK(to, nparts);
      for (size_t k = 0; k < parts.size(); ++k)
        PCU_COMM_PACK(to, parts[k]);
    }
  }
  PCU_Comm_Send();
  /* receiving a global id and associated parts,
     lookup the vertex and classify it on the partition
     model entity for that set of parts */
  while (PCU_Comm_Receive()) {
    Gid gid;
    PCU_COMM_UNPACK(gid);
    int nparts;
    PCU_COMM_UNPACK(nparts);
    Parts residence;
    for (int i = 0; i < nparts; ++i) {
      int part;
      PCU_COMM_UNPACK(part);
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
  int self = PCU_Comm_Self();
  PCU_Comm_Begin();
  APF_ITERATE(GlobalToVert, globalToVert, it) {
    Gid gid = it->first;
    if(gid < 0 )      
        lion_eprint(1, "constructRemotes cgTV4: self=%d,gid=%ld \n",self,gid);
    MeshEntity* vert = it->second;
    Parts residence;
    m->getResidence(vert, residence);
    APF_ITERATE(Parts, residence, rit)
      if (*rit != self) {
        PCU_COMM_PACK(*rit, gid);
        PCU_COMM_PACK(*rit, vert);
      }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    Gid gid;
    PCU_COMM_UNPACK(gid);
    MeshEntity* remote;
    PCU_COMM_UNPACK(remote);
    int from = PCU_Comm_Sender();
    MeshEntity* vert = globalToVert[gid];
    m->addRemote(vert, from, remote);
  }
  // who is not stuck?
  lion_eprint(1, "%d done inside remotes \n",PCU_Comm_Self());
}

void assemble(Mesh2* m, const Gid* conn, int nelem, int etype,
     GlobalToVert& globalToVert)
{
  constructVerts(m, conn, nelem, etype, globalToVert);
  constructElements(m, conn, nelem, etype, globalToVert);
}

void finalise(Mesh2* m, GlobalToVert& globalToVert)
{
  constructResidence(m, globalToVert);
  lion_eprint(1, "%d after residence \n",PCU_Comm_Self());
  constructRemotes(m, globalToVert);
  stitchMesh(m);
  m->acceptChanges();
}
 
void construct(Mesh2* m, const Gid* conn, int nelem, int etype,
    GlobalToVert& globalToVert)
{
  assemble(m, conn, nelem, etype, globalToVert);
  finalise(m, globalToVert);
}

void setCoords(Mesh2* m, const double* coords, int nverts,
    GlobalToVert& globalToVert)
{
  Gid max = getMax(globalToVert);
  Gid total = max + 1;
  int peers = PCU_Comm_Peers();
  int quotient = total / peers;
  int remainder = total % peers;
  int mySize = quotient;
  int self = PCU_Comm_Self();
  if (self == (peers - 1))
    mySize += remainder;
  Gid myOffset = (long)self * quotient;

  /* Force each peer to have exactly mySize verts.
     This means we might need to send and recv some coords */
  double* c = new double[mySize*3];

  Gid start = PCU_Exscan_Long(nverts);

  PCU_Comm_Begin();  // the forced 64 bit math below may not be necessary 
  Gid tmpL=start / quotient; 
  int tmpInt=tmpL;
  int to = std::min(peers - 1, tmpInt);
  tmpL=(to+1)*(long)quotient-start; 
  tmpInt=tmpL;
  int n = std::min(tmpInt, nverts);
  if(n > 1000) {
     Gid sizeToSend=n*3*sizeof(double);
     lion_eprint(1, "setCoords int overflow of: self=%ld,mySize=%ld,total=%ld, n=%ld,to=%ld, quotient=%ld, remainder=%ld start=%ld, peers=%ld, sizeToSend=%ld, nverts=%u \n",self,mySize,total,n,to,quotient,remainder,start,peers,sizeToSend,nverts);
//  Gid peersG = PCU_Comm_Peers();
//    Gid quotientG = total / peersG;
//  Gid remainderG = total % peersG;
//     lion_eprint(1, "setCoords Gid0test: self=%d,mySize=%d,total=%ld, quotientG=%ld, remainderG=%ld,peers=%ld \n",self,mySize,total,quotientG,remainderG,peersG);
}
  PCU_Barrier();

  while (nverts > 0) {
    PCU_COMM_PACK(to, start);
    PCU_COMM_PACK(to, n);
    PCU_Comm_Pack(to, coords, n*3*sizeof(double));

    nverts -= n;
    start += n;
    coords += n*3;
    to = std::min(peers - 1, to + 1);
    n = std::min(quotient, nverts);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    PCU_COMM_UNPACK(start);
    PCU_COMM_UNPACK(n); // |||||| more in-place 64 bit math
    PCU_Comm_Unpack(&c[(start - myOffset) * 3], n*3*sizeof(double));
  }

  /* Tell all the owners of the coords what we need */
  typedef std::vector< std::vector<int> > TmpParts;
  TmpParts tmpParts(mySize);
  PCU_Comm_Begin();
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it) {
    Gid gid = it->first;
    Gid tmpL=gid / quotient;
    int tmpInt=tmpL;
    int to = std::min(peers - 1, tmpInt);
    PCU_COMM_PACK(to, gid);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    Gid gid;
    PCU_COMM_UNPACK(gid);
    Gid from = PCU_Comm_Sender();
    Gid tmpL=gid - myOffset;
    int tmpInt=tmpL;
    tmpParts.at(tmpInt).push_back(from);
  }

  /* Send the coords to everybody who want them */
  PCU_Comm_Begin();
  for (int i = 0; i < mySize; ++i) {
    std::vector<int>& parts = tmpParts[i];
    for (size_t j = 0; j < parts.size(); ++j) {
      int to = parts[j];
      Gid gid = i + myOffset;
      PCU_COMM_PACK(to, gid);
      PCU_Comm_Pack(to, &c[i*3], 3*sizeof(double));
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    Gid gid;
    PCU_COMM_UNPACK(gid);
    double v[3];
    PCU_Comm_Unpack(v, sizeof(v));
    Vector3 vv(v);
    m->setPoint(globalToVert[gid], 0, vv);
  }

  delete [] c;
}

void setMatches(Mesh2* m, const Gid* matches, int nverts,
    GlobalToVert& globalToVert)
{
  Gid max = getMax(globalToVert);
  Gid total = max + 1;
  int peers = PCU_Comm_Peers();
  int quotient = total / peers;
  int remainder = total % peers;
  int mySize = quotient;
  int self = PCU_Comm_Self();
  if (self == (peers - 1))
    mySize += remainder;
  Gid myOffset = (long)self * quotient;

  /* Force each peer to have exactly mySize verts.
     This means we might need to send and recv some matches */
  Gid* c = new Gid[mySize];
  Gid start = PCU_Exscan_Long(nverts);

  PCU_Comm_Begin();

  Gid tmpL=start / quotient; 
  int tmpInt=tmpL;
  int to = std::min(peers - 1, tmpInt);
  tmpL=(to+1)*(long)quotient-start; 
  tmpInt=tmpL;
  int n = std::min(tmpInt, nverts);
  while (nverts > 0) {
    PCU_COMM_PACK(to, start);
    PCU_COMM_PACK(to, n);
    PCU_Comm_Pack(to, matches, n*sizeof(Gid));

    nverts -= n;
    start += n;
    matches += n;
    to = std::min(peers - 1, to + 1);
    n = std::min(quotient, nverts);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    PCU_COMM_UNPACK(start);
    PCU_COMM_UNPACK(n); /// in-place 64
    PCU_Comm_Unpack(&c[(start - myOffset)], n*sizeof(Gid));
  }


  /* Tell all the owners of the matches what we need */
  typedef std::vector< std::vector<int> > TmpParts;
  TmpParts tmpParts(mySize);
  PCU_Comm_Begin();
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it) {
    Gid gid = it->first;  
    int tmpI=gid / quotient;
    int to = std::min(peers - 1,tmpI); 
    PCU_COMM_PACK(to, gid);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    Gid gid;
    PCU_COMM_UNPACK(gid);
    int from = PCU_Comm_Sender();
    tmpParts.at(gid - myOffset).push_back(from);
  }
 
  MeshTag* matchGidTag = m->createLongTag("matchGids", 1);
  /* Send the matches to everybody who wants them */
  PCU_Comm_Begin();
  for (int i = 0; i < mySize; ++i) {
    std::vector<int>& parts = tmpParts[i];
    for (size_t j = 0; j < parts.size(); ++j) {
      int to = parts[j];
      Gid gid = i + myOffset;
      Gid matchGid = c[i];
      PCU_COMM_PACK(to, gid);
      PCU_COMM_PACK(to, matchGid);
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    Gid gid;
    PCU_COMM_UNPACK(gid);
    Gid match;
    PCU_COMM_UNPACK(match);
    PCU_ALWAYS_ASSERT(gid != match);
    PCU_ALWAYS_ASSERT(globalToVert.count(gid));
    m->setLongTag(globalToVert[gid], matchGidTag, &match);
  }

  /* Use the 1D partitioning of global ids to distribute the 
   * entity pointers and their owning ranks. Process 0 will hold
   * the entity pointers and owners for mesh vertex gid [0..quotient),
   * process 1 holds gids [quotient..2*quotient), ...
   */
  PCU_Comm_Begin();
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it) {
    MeshEntity* e = it->second;  // KEJ does not follow this
    Gid gid = it->first;
    int tmpI=gid / quotient;
    int to = std::min(peers - 1,tmpI); 
    PCU_COMM_PACK(to, gid);
    PCU_COMM_PACK(to, e);
  }
  PCU_Comm_Send();
  typedef std::pair< int, apf::MeshEntity* > EntOwnerPtrs;
  typedef std::map< Gid, std::vector< EntOwnerPtrs > > GidPtrs;
  GidPtrs gidPtrs;
  while (PCU_Comm_Receive()) {
    Gid gid;
    PCU_COMM_UNPACK(gid);
    MeshEntity* vert;
    PCU_COMM_UNPACK(vert);
    int owner = PCU_Comm_Sender();
    gidPtrs[gid-myOffset].push_back(EntOwnerPtrs(owner,vert));
  }

  /* Tell the brokers of the matches we need */
  typedef std::pair<Gid,Gid> MatchingPair;
  typedef std::map< MatchingPair, std::vector<Gid> > MatchMap;
  MatchMap matchParts;
  PCU_Comm_Begin();
  APF_CONST_ITERATE(GlobalToVert, globalToVert, it) { //loop over local verts
    Gid gid = it->first;
    Gid matchGid;
    m->getLongTag(it->second, matchGidTag, &matchGid); //get the matched ent gid
    if( matchGid != -1 ) {  // marker for an unmatched vertex
      int tmpI=matchGid / quotient;
      int to = std::min(peers - 1,tmpI);  //broker
      PCU_COMM_PACK(to, gid); // send the local vert gid
      PCU_COMM_PACK(to, matchGid); // and the match gid needed
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    Gid gid;
    PCU_COMM_UNPACK(gid); // request from entity gid
    Gid matchGid;
    PCU_COMM_UNPACK(matchGid); // requesting matched entity gid 
    MatchingPair mp(gid,matchGid);
    int from = PCU_Comm_Sender();
    matchParts[mp].push_back(from); // store a list of the proceses that need the pair (entity gid, match gid)
  }

  /* Send the match pointer and owner process to everybody
   * who wants them */
  PCU_Comm_Begin();
  APF_ITERATE(MatchMap,matchParts,it) {
    MatchingPair mp = it->first;
    Gid gid = mp.first;
    Gid matchGid = mp.second;
    std::vector<Gid> parts = it->second;
    for(size_t i=0; i<parts.size(); i++) {
      const int to = parts[i];
      PCU_COMM_PACK(to, gid);
      PCU_COMM_PACK(to, matchGid);
      size_t numMatches = gidPtrs[matchGid-myOffset].size();
      PCU_COMM_PACK(to, numMatches);
      for( size_t i=0; i<gidPtrs[matchGid-myOffset].size(); i++) {
        EntOwnerPtrs eop = gidPtrs[matchGid-myOffset][i];
        int owner = eop.first;
        PCU_COMM_PACK(to, owner);
        apf::MeshEntity* ent = eop.second;
        PCU_COMM_PACK(to, ent);
      }
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    Gid gid;
    PCU_COMM_UNPACK(gid);
    Gid matchGid;
    PCU_COMM_UNPACK(matchGid);
    size_t numMatches;
    PCU_COMM_UNPACK(numMatches);
    for(size_t i=0; i<numMatches; i++) {
      int owner;
      PCU_COMM_UNPACK(owner);
      MeshEntity* match;
      PCU_COMM_UNPACK(match);
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
