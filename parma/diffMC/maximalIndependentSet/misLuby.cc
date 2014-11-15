#include <mpi.h>
#include <algorithm>
#include <limits>
#include <iostream>
#include <sstream>
#include <iterator>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <locale>

#include <map>
#include <list>
#include <set>

#include "mis.h"

using std::find;
using std::copy;
using std::back_inserter;
using std::map;
using std::list;
using std::make_pair;
using std::min_element;
using std::stable_sort;
using std::unique;
using std::set_intersection;

using std::cout;
using std::vector;
using std::set;

using namespace misLuby;

namespace {
  void debugPrint(const char* msg, bool dbgOverride = false){
    if (dbgOverride) {
      printf("MIS_DEBUG %s", msg);
      fflush(stdout);
    }
  }
}

void misLuby::PartInfo::print() {
  std::string header="adjPartIds ";
  std::stringstream out;
  Print <int, vector<int>::iterator > (out, header, adjPartIds.begin(), 
      adjPartIds.end(), ", ", true);
  header = "net ";
  Print <int, vector<int>::iterator > (out, header, net.begin(), 
      net.end(), ", ", true);
  header = "netAdjParts ";
  Print(out, header, netAdjParts.begin(), netAdjParts.end(), ", ", true);
  out << "randNum " << randNum << "\n";
  std::string msg = out.str();
  PCU_Debug_Print("%s", msg.c_str());
}

void misLuby::Print(std::ostream& ostr, char* dbgMsg,
    vector<adjPart>::iterator itbegin, vector<adjPart>::iterator itend,
    const std::string& delimiter, bool dbg) {
  if (dbg) {
    ostr << "DEBUG " << dbgMsg;
    vector<adjPart>::iterator itr;
    for (itr = itbegin; itr != itend; itr++) {
      ostr << "(" << itr->partId << ", {";
      for (size_t i = 0; i < itr->net.size(); i++) {
        ostr << itr->net[i] << ";";
      }
      ostr << "} )" << delimiter;
    }
    ostr << "\n";
  }
}

void misLuby::Print(std::ostream& ostr, std::string dbgMsg,
    map<int, int>::iterator itbegin, map<int, int>::iterator itend,
    const std::string& delimiter, bool dbg) {
  if (dbg) {
    ostr << "DEBUG " << dbgMsg;
    map<int, int>::iterator itr;
    for (itr = itbegin; itr != itend; itr++) {
      ostr << "(" << itr->first << ";" << itr->second << ")" << delimiter;
    }
    ostr << "\n";
  }
}

void seedRandomNumberGenerator(int seed, int) {
  char msg[256];
  sprintf(msg, "Seeding Random Number Generator with %d\n", seed);
  debugPrint(msg);
  srand(seed);
}

/**
 * @brief generate random number
 * @return 0 on success, non-zero otherwise
 */
int generateRandomNumber() {
  return rand();
}

size_t createNeighborAndGetBufSize(const int destRank) {
  int empty;
  PCU_Comm_Pack(destRank, &empty, 0); // make sure the peer exists
  size_t buffSizeInitial;
  PCU_Comm_Packed(destRank, &buffSizeInitial);
  return buffSizeInitial;
}

int sendAdjNetsToNeighbors(partInfo& part,
    vector<adjPart>& nbNet) {
  PCU_Comm_Begin();

  MIS_ITERATE(vector<int>, part.adjPartIds, adjPartIdItr) {
    const int destRank = *adjPartIdItr;
    MIS_ITERATE(vector<adjPart>, nbNet, nbItr) {
      if (nbItr->partId != *adjPartIdItr) {
        const size_t buffSizeInitial =
          createNeighborAndGetBufSize(destRank);

        //number of bytes in message
        int numIntsInMsg = 5 + nbItr->net.size();
        PCU_COMM_PACK(destRank,numIntsInMsg);

        //pack source part Id
        PCU_COMM_PACK(destRank, part.id);

        //pack destination part Id
        PCU_COMM_PACK( destRank, *adjPartIdItr);

        //adjacent part Id
        PCU_COMM_PACK( destRank, nbItr->partId);

        //adjacent part's random number
        PCU_COMM_PACK( destRank, nbItr->randNum);

        //pack int array
        for (vector<int>::iterator pItr = nbItr->net.begin();
            pItr != nbItr->net.end();
            pItr++) {
          PCU_COMM_PACK( destRank, *pItr);
        }

        //sanity check
        size_t buffSizeFinal;
        PCU_Comm_Packed(destRank, &buffSizeFinal);
        const size_t numIntsPacked =
          (buffSizeFinal - buffSizeInitial) / sizeof (int);
        MIS_FAIL_IF((size_t)numIntsInMsg != numIntsPacked,
            "number of ints packed does not match msg header");
      }
    }
  }
  return 0;
}



int sendIntsToNeighbors(vector<partInfo>& parts,
    vector<int>*& msg, set<int>& destRanks, int tag) {
  const int rank = PCU_Comm_Self();
  PCU_Comm_Start(PCU_GLOBAL_METHOD);

  int partIdx = 0;
  MIS_ITERATE(vector<partInfo>, parts, partItr) {
    MIS_ITERATE(vector<int>, partItr->adjPartIds, adjPartIdItr) {
      const int destRank = *adjPartIdItr;
      destRanks.insert(destRank);

      const size_t buffSizeInitial = createNeighborAndGetBufSize(destRank);

      //number of bytes in message
      int numIntsInMsg = 4 + msg[partIdx].size();
      PCU_COMM_PACK( destRank, numIntsInMsg);

      //pack msg tag
      PCU_COMM_PACK( destRank, tag);

      //pack source part Id
      PCU_COMM_PACK( destRank, partItr->id);

      //pack destination part Id
      PCU_COMM_PACK( destRank, *adjPartIdItr);

      //pack int array
      if (msg[partIdx].size() == 0) {
        char dbgMsg[256];
        sprintf(dbgMsg,
            "[%d] (%d) %s packing empty integer array\n",
            rank, partItr->id, __FUNCTION__);
        debugPrint(dbgMsg);
      }
      for ( vector<int>::iterator pItr = msg[partIdx].begin();
          pItr != msg[partIdx].end();
          pItr++ ) {
        PCU_COMM_PACK( destRank, *pItr);
      }

      //sanity check
      size_t buffSizeFinal;
      PCU_Comm_Packed(destRank, &buffSizeFinal);
      const size_t numIntsPacked =
        (buffSizeFinal - buffSizeInitial) / sizeof (int);
      MIS_FAIL_IF((size_t)numIntsInMsg != numIntsPacked,
          "number of ints packed does not match msg header");
    }
    partIdx++;
  }
  return 0;
}

void unpackInts(vector<partInfo>& parts, vector<int>*& msg,
    const size_t numIntsInMsg, int tag) {
  const int rank = PCU_Comm_Self();
  int srcRank;
  PCU_Comm_From(&srcRank);

  size_t numIntsUnpacked = 0;
  //unpack msg tag
  int inTag;
  PCU_COMM_UNPACK(inTag);
  MIS_FAIL_IF(tag != inTag, "tags do not match");
  numIntsUnpacked++;

  //unpack source part Id
  int srcPartId;
  PCU_COMM_UNPACK(srcPartId);
  assert(srcRank == srcPartId);
  numIntsUnpacked++;

  //unpack destination part Id
  int destPartId;
  PCU_COMM_UNPACK(destPartId);
  assert(rank == destPartId);
  numIntsUnpacked++;

  //find index to store message
  size_t partIdx = 0;
  for (; partIdx < parts.size(); partIdx++) {
    if (parts[partIdx].id == destPartId) {
      break;
    }
  }
  assert(partIdx < parts.size());

  //unpack int array
  const size_t buffSize = (numIntsInMsg - numIntsUnpacked) * sizeof (int);
  if (buffSize == 0) {
    char dbgMsg[256];
    sprintf(dbgMsg, "[%d] (%d) %s unpacking empty integer array\n",
        rank, destPartId, __FUNCTION__);
    debugPrint(dbgMsg);
  }

  int buff;
  for(size_t i=0; i < numIntsInMsg - numIntsUnpacked; i++) {
    PCU_COMM_UNPACK(buff);
    msg[partIdx].push_back(buff);
  }
}

void recvIntsFromNeighbors(vector<partInfo>& parts,
    vector<int>*& msg, set<int>& srcRanks, int, int tag) {
  while(PCU_Comm_Listen()) {
    size_t msgSz;
    PCU_Comm_Received(&msgSz);
    const size_t numIntsInBuff =  msgSz / sizeof (int);
    size_t numIntsProcessed = 0;

    int srcRank;
    PCU_Comm_From(&srcRank);
    srcRanks.insert(srcRank);

    do {
      int numIntsPacked;
      PCU_COMM_UNPACK(numIntsPacked);
      numIntsProcessed += numIntsPacked;
      unpackInts(parts, msg, numIntsPacked - 1, tag);
    } while (numIntsProcessed != numIntsInBuff);
  }
}

int sendNetToNeighbors(partInfo& part) {
  PCU_Comm_Begin();
  MIS_ITERATE(vector<int>, part.adjPartIds, adjPartIdItr) {
    const int destRank = *adjPartIdItr;
    PCU_COMM_PACK(destRank, part.id); //FIXME redundant
    PCU_COMM_PACK(destRank, *adjPartIdItr); //FIXME redundant
    PCU_COMM_PACK(destRank, part.randNum);
    size_t netSz = part.net.size();
    PCU_COMM_PACK(destRank, netSz);
    MIS_ITERATE(vector<int>, part.net, pItr)
      PCU_COMM_PACK(destRank, *pItr);
  }
  return 0;
}

void unpackNet(partInfo& part, vector<adjPart>& msg) {
  int srcPartId;
  PCU_COMM_UNPACK(srcPartId);
  assert(PCU_Comm_Sender() == srcPartId);

  int destPartId;
  PCU_COMM_UNPACK(destPartId);
  assert(PCU_Comm_Self() == destPartId);

  int randNum;
  PCU_COMM_UNPACK(randNum);

  adjPart ap;
  ap.partId = srcPartId;
  ap.randNum = randNum;

  size_t arrayLen;
  PCU_COMM_UNPACK(arrayLen);
  assert(arrayLen > 0);

  int buff;
  for(size_t i=0; i < arrayLen; i++) {
    PCU_COMM_UNPACK(buff);
    ap.net.push_back(buff);
  }

  msg.push_back(ap);
}

void recvNetsFromNeighbors(partInfo& part, vector<adjPart>& msg) {
  while( PCU_Comm_Listen() )
    unpackNet(part, msg);
}

void unpackAdjPart(vector<partInfo>& parts,
    vector<adjPart>*& msg, const size_t numIntsInMsg) {
  const int rank = PCU_Comm_Self();
  int srcRank;
  PCU_Comm_From(&srcRank);

  size_t numIntsUnpacked = 0;
  //unpack source part Id
  int srcPartId;
  PCU_COMM_UNPACK(srcPartId);
  assert(srcRank == srcPartId);
  numIntsUnpacked++;

  //unpack destination part Id
  int destPartId;
  PCU_COMM_UNPACK(destPartId);
  assert(rank == destPartId);
  numIntsUnpacked++;

  //unpack adjacent part Id
  int adjPartId;
  PCU_COMM_UNPACK(adjPartId);
  numIntsUnpacked++;

  //unpack random number
  int randNum;
  PCU_COMM_UNPACK(randNum);
  numIntsUnpacked++;

  //find index to store message
  size_t partIdx = 0;
  for (; partIdx < parts.size(); partIdx++) {
    if (parts[partIdx].id == destPartId) {
      break;
    }
  }
  assert(partIdx < parts.size());

  adjPart ap;
  ap.partId = adjPartId;
  ap.randNum = randNum;

  //unpack int array
  const size_t buffSize = (numIntsInMsg - numIntsUnpacked) * sizeof (int);
  assert(buffSize > 0);

  int buff;
  for(size_t i=0; i < numIntsInMsg - numIntsUnpacked; i++) {
    PCU_COMM_UNPACK(buff);
    ap.net.push_back(buff);
  }

  msg[partIdx].push_back(ap);
}

void recvAdjNetsFromNeighbors(vector<partInfo>& parts,
    vector<adjPart>*& msg) {

  while( PCU_Comm_Listen() ) {
    size_t msgSz;
    PCU_Comm_Received(&msgSz);
    const size_t numIntsInBuff = msgSz / sizeof (int);
    size_t numIntsProcessed = 0;

    do {
      int numIntsPacked;
      PCU_COMM_UNPACK(numIntsPacked);
      numIntsProcessed += numIntsPacked;
      unpackAdjPart(parts, msg, numIntsPacked - 1);
    } while (numIntsProcessed != numIntsInBuff);
  }
}

int minRandNum(netAdjItr first, netAdjItr last) {
  int min = INT_MAX;
  for (netAdjItr itr = first; itr != last; itr++) {
    if (itr->second < min) {
      min = itr->second;
    }
  }
  return min;
}

void getNetAdjPartIds(netAdjItr first, netAdjItr last, vector<int>& ids) {
  for (netAdjItr itr = first; itr != last; itr++) {
    ids.push_back(itr->first);
  }
}

void setNodeStateInGraph(partInfo& part) {
  part.isInNetGraph = true;
  part.isInMIS = false;
}

/**
 * @brief add parts as neighbors in the net-graph
 * @remark Intersect the local net with each of the neighbors nets.  If the
 *         intersection is not empty then add the neighbor to the netAdjPart
 *         list.
 * @param nbNet (In) neighbor parts and their nets
 */
void partInfo::addNetNeighbors(vector<adjPart>& nbNet) {
  stable_sort(this->net.begin(), this->net.end());

  //loop over each neighbor to the part
  MIS_ITERATE(vector<adjPart>, nbNet, nbNetItr) {
    stable_sort(nbNetItr->net.begin(), nbNetItr->net.end());

    vector<int> netIntersection;
    set_intersection(nbNetItr->net.begin(), nbNetItr->net.end(),
        this->net.begin(), this->net.end(),
        back_inserter(netIntersection));

    if (netIntersection.size() > 0) {
      this->netAdjParts[nbNetItr->partId] = nbNetItr->randNum;
    }
  }
}

/**
 * @brief update the neighbor part list after changing the list of
 *        net-graph neighbors
 */
void partInfo::updateNeighbors(){

  vector<int>::iterator findItr;
  MIS_ITERATE(mapIntInt, netAdjParts, netAdjPartItr) {
    findItr = find(adjPartIds.begin(),
                   adjPartIds.end(),
                   netAdjPartItr->first);
    if (findItr == adjPartIds.end())
      adjPartIds.push_back(netAdjPartItr->first);
  }
  stable_sort(adjPartIds.begin(), adjPartIds.end());
}

/**
 * @brief for each part determine with other parts are net-graph neighbors
 * @remark 2nd-adjacent net neighbors are found with the second round of
 *         communications
 *
 * The netgraph is constructed with one node for each part and one edge between
 * two parts if the intersection of their nets is not empty.  For example,
 * consider three parts:
 *
 *  a  b  c
 *
 * where a's net is b, b's net is a, and c's net is a.  In the figure below the
 * net for each node is noted with dashes:
 *
 *  a b c
 * a---
 * b---
 * c  ---
 *
 * The netgraph is then:
 *  a-b-c
 *  \___/
 *
 * With nets defined as N,S,E,W neighbors each node in the netgraph has an edge
 * to its W,N,E,S neighbors through a first adjacency; i.e. an edge between
 * node u and v is created as node u has v in it's net and visa versa.  Each
 * node also has an edge to it's NW, NE, SE, NE, WW, NN, EE, and SS neighbors
 * through second adjacenies; i.e. node u and v both have node w!=u,v, in
 * their nets.
 *
 * The figure below depicts these adjacencies with a 0 indicating the given
 * node, a 1 indicating a first-adjacency, and a 2 indicating a second
 * adjacency.
 *
 *     2
 *   2 1 2
 * 2 1 0 1 2
 *   2 1 2
 *     2
 *
 * The second-adjacent neighbors are defined by two parts that don't contain
 * each other in their nets, but do have intersecting nets.
 *
 * @param parts (InOut) vector of local parts
 * @return 0 on success, non-zero otherwise
 */
int constructNetGraph(partInfo& part) {

  int rank;
  PCU_Comm_Rank(&rank);

  stable_sort(part.adjPartIds.begin(), part.adjPartIds.end());
  setNodeStateInGraph(part);

  // send net to each neighbor
  int ierr = sendNetToNeighbors(part);
  if (ierr != 0) return ierr;
  PCU_Comm_Send();
  vector<adjPart> nbNet;
  recvNetsFromNeighbors(part, nbNet);

  // get the random num associated with adjacent nets
  part.addNetNeighbors(nbNet);

  ierr = sendAdjNetsToNeighbors(part, nbNet);
  if (ierr != 0) return ierr;
  PCU_Comm_Send();
  vector<adjPart>* nbAdjNets = new vector<adjPart>[1]; //FIXME
  vector<PartInfo> parts; parts.push_back(part); //FIXME
  recvAdjNetsFromNeighbors(parts, nbAdjNets); //FIXME

  // get the random num associated with adjacent nets
  part.addNetNeighbors(nbAdjNets[0]); //FIXME
  delete [] nbAdjNets;

  //add new neighbors based on changes to net-graph neighbors
  MIS_ITERATE(vector<partInfo>, parts, partItr)
    partItr->updateNeighbors();

  return 0;
}

void setRandomNum(partInfo& part) {
  part.randNum = generateRandomNumber();
  // don't select self nets until all other nets are selected
  if ( 1 == part.net.size() )
    part.randNum = std::numeric_limits<int>::max();
}

void mis_init(int randNumSeed, int debugMode, const char* maj,
    const char* min, const char* patch) {
  char msg[256];
  sprintf(msg, "MIS version %s.%s.%s\n", maj, min, patch);
  debugPrint(msg);
  seedRandomNumberGenerator(randNumSeed, debugMode);
}

void misFinalize() {
}

bool doSetsMatch(set<int> a, set<int> b) {
  set<int>::iterator aItr = a.begin();
  set<int>::iterator bItr = b.begin();
  while (aItr != a.end() && bItr != b.end()) {
    if(*aItr != *bItr) return false;
    aItr++;
    bItr++;
  }
  return true;
}


int mis(partInfo& part, vector<int>& mis, bool randNumsPredefined) {
  assert(PCU_Comm_Initialized());

  if (false == randNumsPredefined)
    setRandomNum(part);

  int ierr = constructNetGraph(part);
  if (ierr != 0) return ierr;

  vector<partInfo> parts; parts.push_back(part); //FIXME

  const int numParts = 1;
  vector<int>* nodesRemoved = new vector<int>[numParts];
  vector<int>* nodesToRemove = new vector<int>[numParts];
  vector<int>* rmtNodesToRemove = new vector<int>[numParts];

  int totalNodesAdded = 0;
  int loopCount = 0;
  int tag = 0;
  do {
    int numNodesAdded = 0;
    const int minRand = 
      minRandNum(part.netAdjParts.begin(), part.netAdjParts.end());
    if (true == part.isInNetGraph &&
        false == part.isInMIS &&
        part.randNum < minRand) {
      part.isInMIS = true;
      mis.push_back(part.id);
      ++numNodesAdded;
    }

    const int partIdx = 0;
    if (true == part.isInMIS) {
      nodesToRemove[partIdx].reserve(part.netAdjParts.size() + 1);
      getNetAdjPartIds(part.netAdjParts.begin(),
          part.netAdjParts.end(), nodesToRemove[partIdx]);
      nodesToRemove[partIdx].push_back(part.id);
    }

    set<int> destRanks;
    tag++;
    ierr = sendIntsToNeighbors(parts, nodesToRemove, destRanks, tag);
    MIS_FAIL_IF(ierr != 0, "sendIntsToNeighbors failed");
    PCU_Comm_Send();

    set<int> srcRanks;
    recvIntsFromNeighbors(parts, rmtNodesToRemove, srcRanks, loopCount, tag);
    MIS_FAIL_IF(destRanks.size() != srcRanks.size(),
        "recvIntsFromNeighbors failed");
    MIS_FAIL_IF(!doSetsMatch(srcRanks, destRanks), "ranks do not match");

    if (true == part.isInNetGraph) {
      if (true == part.isInMIS ||
          find(rmtNodesToRemove[partIdx].begin(),
            rmtNodesToRemove[partIdx].end(),
            part.id) != rmtNodesToRemove[partIdx].end()) {
        part.netAdjParts.clear();
        part.isInNetGraph = false;
        nodesRemoved[partIdx].push_back(part.id);
      } else {
        MIS_ITERATE(vector<int>, rmtNodesToRemove[partIdx], rmvNodeItr)
          part.netAdjParts.erase(*rmvNodeItr);
      }
    }
    nodesToRemove[partIdx].clear();
    rmtNodesToRemove[partIdx].clear();

    destRanks.clear();
    tag++;
    ierr = sendIntsToNeighbors(parts, nodesRemoved, destRanks, tag);
    MIS_FAIL_IF(ierr != 0, "sendIntsToNeighbors failed");
    PCU_Comm_Send();

    srcRanks.clear();
    recvIntsFromNeighbors(parts, rmtNodesToRemove, srcRanks, loopCount, tag);
    MIS_FAIL_IF(destRanks.size() != srcRanks.size(),
        "recvIntsFromNeighbors failed");
    MIS_FAIL_IF(!doSetsMatch(srcRanks, destRanks), "ranks do not match");

    if (true == part.isInNetGraph)
      MIS_ITERATE(vector<int>, rmtNodesToRemove[partIdx], rmvNodeItr)
        part.netAdjParts.erase(*rmvNodeItr);
    nodesRemoved[partIdx].clear();
    rmtNodesToRemove[partIdx].clear();

    totalNodesAdded+=numNodesAdded;
    PCU_Add_Ints(&totalNodesAdded, 1);
    loopCount++;
  } while (totalNodesAdded > 0);

  char dbgMsg[512];
  sprintf(dbgMsg, "Number of mis loops: %d\n", loopCount);
  debugPrint(dbgMsg);

  delete [] nodesRemoved;
  delete [] nodesToRemove;
  delete [] rmtNodesToRemove;

  return 0;
}

void sendMisStatusToNeighbors(vector<partInfo>& parts) {

  PCU_Comm_Start(PCU_GLOBAL_METHOD);

  //pack
  int partIdx = 0;
  MIS_ITERATE(vector<partInfo>, parts, partItr) {
    MIS_ITERATE(vector<int>, partItr->adjPartIds, adjPartIdItr) {
      const int destRank = *adjPartIdItr;

      createNeighborAndGetBufSize(destRank);

      //number of integers in message
      int numIntsInMsg = 3; //TODO remove this
      PCU_COMM_PACK( destRank, numIntsInMsg); //TODO remove this

      //pack destination part Id
      PCU_COMM_PACK( destRank, *adjPartIdItr);

      //pack source part Id or -1 if source part is not in MIS
      int srcPartId = -1;
      vector<int>::iterator it;
      it = find(partItr->net.begin(), partItr->net.end(), *adjPartIdItr);
      if (it != partItr->net.end() && true == partItr->isInMIS) {
        srcPartId = partItr->id;
      } else {
        srcPartId = -1;
      }
      PCU_COMM_PACK( destRank, srcPartId);
    }
    partIdx++;
  }
}

void recvMisStatusFromNeighbors(vector<partInfo>& parts,
    vector<int>& neighborsInMIS) {

  const int rank = PCU_Comm_Self();

  neighborsInMIS.resize(parts.size(), -1);

  assert(neighborsInMIS.size() == parts.size());
  MIS_ITERATE(vector<int>, neighborsInMIS, vItr)
    assert(-1 == *vItr);

  while (PCU_Comm_Listen()) {
    size_t msgSz;
    PCU_Comm_Received(&msgSz);
    const size_t numIntsInBuff = msgSz / sizeof (int);
    size_t numIntsProcessed = 0;

    int srcRank;
    PCU_Comm_From(&srcRank);
    do {
      size_t numIntsPacked;
      size_t numIntsUnpacked = 0;
      PCU_COMM_UNPACK(numIntsPacked);
      numIntsUnpacked++;
      numIntsProcessed += numIntsPacked;

      //unpack destination part Id
      int destPartId;
      PCU_COMM_UNPACK(destPartId);
      assert(rank == destPartId);
      numIntsUnpacked++;

      //find index to store message
      size_t partIdx = 0;
      for (; partIdx < parts.size(); partIdx++) {
        if (parts[partIdx].id == destPartId) {
          break;
        }
      }
      assert(partIdx < parts.size());

      //unpack source part Id
      int srcPartId;
      PCU_COMM_UNPACK(srcPartId);
      if (-1 != srcPartId) {
        assert(srcRank == srcPartId);
        const int val = neighborsInMIS.at(partIdx);
        char msg[512];
        sprintf(msg,"[%d] (%d) mergeTarget=%d\n", rank,
            parts[partIdx].id, srcPartId);
        debugPrint(msg);
        if (-1 != neighborsInMIS.at(partIdx)) {
          if (neighborsInMIS.at(partIdx) != srcPartId) {
            printf("[%d] srcPartId=%d does not match neighborsInMIS.at(%ld) = %d\n",
                rank, srcPartId, (long)partIdx, val);
          }
        }
        assert(-1 == neighborsInMIS.at(partIdx));
        neighborsInMIS[partIdx] = srcPartId;
      }
      numIntsUnpacked++;

      assert(numIntsPacked == numIntsUnpacked);
    } while (numIntsProcessed != numIntsInBuff);
  }
}

void getMergeTargets(const int, const int,
    vector<partInfo>& parts, vector<int>&, map<int, int>& mergeTargets) {
  if(!PCU_Comm_Initialized()) PCU_Comm_Init();

  sendMisStatusToNeighbors(parts);
  PCU_Comm_Send();
  vector<int> mTgts;
  recvMisStatusFromNeighbors(parts, mTgts);

  int partIdx = 0;
  MIS_ITERATE(vector<misLuby::partInfo>, parts, partItr) {
    if (-1 != mTgts.at(partIdx))
      mergeTargets[partItr->id] = mTgts.at(partIdx);
    partIdx++;
  }
}
