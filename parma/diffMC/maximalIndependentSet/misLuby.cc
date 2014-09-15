#include <algorithm>
#include <limits>
#include <iostream>
#include <iterator>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>
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

void misLuby::Print(std::ostream& ostr, char* dbgMsg, 
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



inline int GetGlobalPartId(int rank, int partIdx, int numParts) {
  return rank * numParts + partIdx;
}

inline int GetOwningProcessRank(int globalPartId, int numParts) {
  assert(numParts >= 1);
  return globalPartId / numParts;
}

void seedRandomNumberGenerator(int seed, int) {
  char msg[256];
  sprintf(msg, "Seeding Random Number Generator with %d\n", seed);
  debugPrint(msg);
  srand(seed);
}

/**
 * @brief generate n random numbers
 * @param randNums (InOut) pre-allocated vector of size n
 * @return 0 on success, non-zero otherwise
 */
int generateRandomNumbers(vector<int>& randNums) {//, bool seedGenerator) {
  for (size_t i = 0; i < randNums.size(); i++) {
    randNums[i] = rand();
  }
  return 0;
}

size_t createNeighborAndGetBufSize(const int destRank) {
  int empty;
  PCU_Comm_Pack(destRank, &empty, 0); // make sure the peer exists
  size_t buffSizeInitial;
  PCU_Comm_Packed(destRank, &buffSizeInitial);
  return buffSizeInitial;
}

int sendAdjNetsToNeighbors(vector<partInfo>& parts, 
    vector<adjPart>*& nbNet) {  

  PCU_Comm_Start(PCU_GLOBAL_METHOD);

  int partIdx = 0;
  MIS_ITERATE(vector<partInfo>, parts, partItr) {
    MIS_ITERATE(vector<int>, partItr->adjPartIds, adjPartIdItr) {
      const int destRank = GetOwningProcessRank(*adjPartIdItr, parts.size());
      MIS_ITERATE(vector<adjPart>, nbNet[partIdx], nbItr) {
        if (nbItr->partId != *adjPartIdItr) {
          const size_t buffSizeInitial = 
            createNeighborAndGetBufSize(destRank);

          //number of bytes in message
          int numIntsInMsg = 5 + nbItr->net.size();
          PCU_COMM_PACK(destRank,numIntsInMsg);

          //pack source part Id
          PCU_COMM_PACK(destRank, partItr->id);

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
    partIdx++;
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
      const int destRank = GetOwningProcessRank(*adjPartIdItr, parts.size());
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
  assert(srcRank == GetOwningProcessRank(srcPartId, parts.size()));
  numIntsUnpacked++;

  //unpack destination part Id
  int destPartId;
  PCU_COMM_UNPACK(destPartId);
  assert(rank == GetOwningProcessRank(destPartId, parts.size()));
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

int sendNetToNeighbors(vector<partInfo>& parts) {

  PCU_Comm_Start(PCU_GLOBAL_METHOD); //TODO update this

  int partIdx = 0;
  MIS_ITERATE(vector<partInfo>, parts, partItr) {
    MIS_ITERATE(vector<int>, partItr->adjPartIds, adjPartIdItr) {
      const int destRank = GetOwningProcessRank(*adjPartIdItr, parts.size());
      //pack source part Id
      PCU_COMM_PACK( destRank, partItr->id);
      //pack destination part Id
      PCU_COMM_PACK( destRank, *adjPartIdItr);
      //part's random number
      PCU_COMM_PACK( destRank, partItr->randNum);
      //pack size of int array and its contents
      size_t netSz = partItr->net.size();
      PCU_COMM_PACK( destRank, netSz);
      MIS_ITERATE(vector<int>, partItr->net, pItr)
        PCU_COMM_PACK( destRank, *pItr);    
    }
    partIdx++;
  }
  return 0;
}

void unpackNet(vector<partInfo>& parts, vector<adjPart>*& msg) {
  int rank;
  PCU_Comm_Rank(&rank);

  int srcRank;
  PCU_Comm_From(&srcRank);

  //unpack source part Id
  int srcPartId;
  PCU_COMM_UNPACK(srcPartId);
  assert(srcRank == GetOwningProcessRank(srcPartId, parts.size())); 

  //unpack destination part Id
  int destPartId;
  PCU_COMM_UNPACK(destPartId);
  assert(rank == GetOwningProcessRank(destPartId, parts.size()));

  //unpack random number
  int randNum;
  PCU_COMM_UNPACK(randNum);

  //find index to store message
  size_t partIdx = 0;
  for (; partIdx < parts.size(); partIdx++) {
    if (parts[partIdx].id == destPartId) {
      break;
    }
  }
  if (partIdx >= parts.size()) {
    fprintf(stderr, "[%d] %s destPartId= %d srcPartId= %d randNum= %d\n",
        rank, __FUNCTION__, destPartId, srcPartId, randNum);
    exit(1);
  }
  assert(partIdx < parts.size());

  adjPart ap;
  ap.partId = srcPartId;
  ap.randNum = randNum;

  //unpack size of int array
  size_t arrayLen;
  PCU_COMM_UNPACK(arrayLen);
  assert(arrayLen > 0);

  int buff;
  for(size_t i=0; i < arrayLen; i++) {
    PCU_COMM_UNPACK(buff);    
    ap.net.push_back(buff);
  }    

  msg[partIdx].push_back(ap);
}

void recvNetsFromNeighbors(vector<partInfo>& parts, 
    vector<adjPart>*& msg) {
  while( PCU_Comm_Listen() )
    unpackNet(parts, msg);
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
  assert(srcRank == GetOwningProcessRank(srcPartId, parts.size()));
  numIntsUnpacked++;

  //unpack destination part Id
  int destPartId;
  PCU_COMM_UNPACK(destPartId);
  assert(rank == GetOwningProcessRank(destPartId, parts.size()));
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

void setNodeStateInGraph(vector<partInfo>& parts) {
  MIS_ITERATE(vector<partInfo>, parts, partItr) {
    partItr->isInNetGraph = true;
    partItr->isInMIS = false;
  }
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
 * @param parts (InOut) vector of local parts
 * @return 0 on success, non-zero otherwise
 */
int constructNetGraph(vector<partInfo>& parts) {
   
  int rank;
  PCU_Comm_Rank(&rank);    

  MIS_ITERATE(vector<partInfo>, parts, partItr)
    stable_sort(partItr->adjPartIds.begin(), partItr->adjPartIds.end());
  setNodeStateInGraph(parts);

  // send net to each neighbor
  int ierr = sendNetToNeighbors(parts);
  if (ierr != 0) return ierr;
  PCU_Comm_Send();
  vector<adjPart>* nbNet = new vector<adjPart>[parts.size()];
  recvNetsFromNeighbors(parts, nbNet);

  // get the random num associated with adjacent nets
  int partIdx = 0;
  MIS_ITERATE(vector<partInfo>, parts, pItr) 
     pItr->addNetNeighbors(nbNet[partIdx++]);

  ierr = sendAdjNetsToNeighbors(parts, nbNet);
  if (ierr != 0) return ierr;
  PCU_Comm_Send();
  delete [] nbNet;
  vector<adjPart>* nbAdjNets = new vector<adjPart>[parts.size()];
  recvAdjNetsFromNeighbors(parts, nbAdjNets);

  // get the random num associated with adjacent nets
  partIdx = 0;
  MIS_ITERATE(vector<partInfo>, parts, partItr) 
    partItr->addNetNeighbors(nbAdjNets[partIdx++]);
  delete [] nbAdjNets;

  //add new neighbors based on changes to net-graph neighbors
  MIS_ITERATE(vector<partInfo>, parts, partItr)
    partItr->updateNeighbors();

  return 0;
}

void setRandomNums(int, vector<partInfo>& parts) {
  const int numParts = parts.size();
  int* localRandNums = new int[numParts];
  vector<int> randNums(parts.size());
  generateRandomNumbers(randNums);

  int partIdx = 0;
  MIS_ITERATE(vector<partInfo>, parts, partItr) {
    if ( 1 == partItr->net.size() ) {
      // don't select self nets until all other nets are selected
      partItr->randNum = std::numeric_limits<int>::max();
    } else {
      partItr->randNum = localRandNums[partIdx++];
    }
  }
  delete [] localRandNums;
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


int mis(const int rank, const int, vector<partInfo>& parts, 
    vector<int>& mis, bool randNumsPredefined) {
  if(!PCU_Comm_Initialized()) PCU_Comm_Init();
  MPI_Barrier(MPI_COMM_WORLD);

  if (false == randNumsPredefined) 
    setRandomNums(rank, parts);

  int ierr = constructNetGraph(parts);
  if (ierr != 0) return ierr;

  vector<int>* nodesRemoved = new vector<int>[parts.size()];
  vector<int>* nodesToRemove = new vector<int>[parts.size()];
  vector<int>* rmtNodesToRemove = new vector<int>[parts.size()];

  char dbgMsg[512];
  int totalNodesAdded = 0;
  int loopCount = 0;
  int tag = 0;
  do {
    int numNodesAdded = 0;
    MIS_ITERATE(vector<partInfo>, parts, partItr) {
      const int minRand = minRandNum(partItr->netAdjParts.begin(), 
          partItr->netAdjParts.end());
      if (true == partItr->isInNetGraph &&
          false == partItr->isInMIS &&
          partItr->randNum < minRand) {
        partItr->isInMIS = true;
        mis.push_back(partItr->id);
        ++numNodesAdded;
      }
    }

    int partIdx = 0;
    MIS_ITERATE(vector<partInfo>, parts, partItr) {
      if (true == partItr->isInMIS) {
        nodesToRemove[partIdx].reserve(partItr->netAdjParts.size() + 1);
        getNetAdjPartIds(partItr->netAdjParts.begin(), 
            partItr->netAdjParts.end(), nodesToRemove[partIdx]);
        nodesToRemove[partIdx].push_back(partItr->id);
      }
      partIdx++;
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

    partIdx = 0;
    MIS_ITERATE(vector<partInfo>, parts, partItr) {
      if (true == partItr->isInNetGraph) {
        if (true == partItr->isInMIS ||
            find(rmtNodesToRemove[partIdx].begin(),
              rmtNodesToRemove[partIdx].end(),
              partItr->id) != rmtNodesToRemove[partIdx].end()) {
          partItr->netAdjParts.clear();
          partItr->isInNetGraph = false;
          nodesRemoved[partIdx].push_back(partItr->id);
        } else {
          MIS_ITERATE(vector<int>, rmtNodesToRemove[partIdx], rmvNodeItr) 
            partItr->netAdjParts.erase(*rmvNodeItr);
        }
      }
      nodesToRemove[partIdx].clear();
      rmtNodesToRemove[partIdx].clear();
      partIdx++;
    }

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

    partIdx = 0;
    MIS_ITERATE(vector<partInfo>, parts, partItr) {
      if (true == partItr->isInNetGraph) {
        MIS_ITERATE(vector<int>, rmtNodesToRemove[partIdx], rmvNodeItr) 
          partItr->netAdjParts.erase(*rmvNodeItr);
      }
      nodesRemoved[partIdx].clear();
      rmtNodesToRemove[partIdx].clear();
      partIdx++;
    }

    MPI_Allreduce(&numNodesAdded, &totalNodesAdded, 1, MPI_INT, MPI_SUM, 
        MPI_COMM_WORLD);
    loopCount++;

    sprintf(dbgMsg, "loopCount = %d\n", loopCount);
    debugPrint(dbgMsg);
  } while (totalNodesAdded > 0);

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
      const int destRank = GetOwningProcessRank(*adjPartIdItr, parts.size());

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
      assert(rank == GetOwningProcessRank(destPartId, parts.size()));
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
        assert(srcRank == GetOwningProcessRank(srcPartId, parts.size()));
        const int val = neighborsInMIS.at(partIdx);
        char msg[512];
        sprintf(msg,"[%d] (%d) mergeTarget=%d\n", rank, 
            parts[partIdx].id, srcPartId);
        debugPrint(msg);
        if (-1 != neighborsInMIS.at(partIdx)) {
          if (neighborsInMIS.at(partIdx) != srcPartId) {
            printf("[%d] srcPartId=%d does not match neighborsInMIS.at(%lu) = %d\n", 
                rank, srcPartId, partIdx, val);
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
  MPI_Barrier(MPI_COMM_WORLD);

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
