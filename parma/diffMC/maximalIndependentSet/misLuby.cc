#include <mpi.h>
#include <algorithm>
#include <limits>
#include <iostream>
#include <sstream>
#include <iterator>
#include <time.h>
#include <stdlib.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <locale>

#include <map>
#include <list>
#include <set>

#include "mis.h"
#include "mersenne_twister.h"

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
  int generateRandomNumber(pcu::PCU *PCUObj) {
    return (int) mersenne_twister(PCUObj);
  }

  void setRandomNum(partInfo& part, pcu::PCU *PCUObj) {
    part.randNum = generateRandomNumber(PCUObj);
    // don't select self nets until all other nets are selected
    if ( 1 == part.net.size() )
      part.randNum = std::numeric_limits<unsigned>::max();
  }

  void removeNodes(partInfo& p, vector<int>& nodes) {
    if (true == p.isInNetGraph)
      MIS_ITERATE(vector<int>, nodes, nodeItr)
        p.netAdjParts.erase(*nodeItr);
  }

  int sendAdjNetsToNeighbors(partInfo& part, vector<adjPart>& nbNet, pcu::PCU *PCUObj) {
    PCUObj->Begin();

    MIS_ITERATE(vector<int>, part.adjPartIds, adjPartIdItr) {
      const int destRank = *adjPartIdItr;
      MIS_ITERATE(vector<adjPart>, nbNet, nbItr) {
        if (nbItr->partId != *adjPartIdItr) {
          //pack destination part Id
          PCUObj->Pack(destRank, *adjPartIdItr);

          //adjacent part Id
          PCUObj->Pack(destRank, nbItr->partId);

          //adjacent part's random number
          PCUObj->Pack(destRank, nbItr->randNum);

          //net size
          size_t n = nbItr->net.size();
          PCUObj->Pack(destRank, n);

          //net
          for (vector<int>::iterator pItr = nbItr->net.begin();
              pItr != nbItr->net.end();
              pItr++) {
            PCUObj->Pack( destRank, *pItr);
          }
        }
      }
    }
    return 0;
  }

  void sendIntsToNeighbors(partInfo& part, vector<int>& msg, int tag, pcu::PCU *PCUObj) {
    PCUObj->Begin();
    MIS_ITERATE(vector<int>, part.adjPartIds, adjPartIdItr) {
      const int destRank = *adjPartIdItr;

      //pack msg tag
      PCUObj->Pack(destRank, tag);

      //pack destination part Id
      PCUObj->Pack(destRank, *adjPartIdItr);

      //pack array length
      size_t n = msg.size();
      PCUObj->Pack(destRank, n);

      //pack int array
      for ( vector<int>::iterator pItr = msg.begin();
          pItr != msg.end();
          pItr++ ) {
        PCUObj->Pack(destRank, *pItr);
      }
    }
  }

  void unpackInts(vector<int>& msg, int tag, pcu::PCU *PCUObj) {
    const int rank = PCUObj->Self();

    //unpack msg tag
    int inTag;
    PCUObj->Unpack(inTag);
    MIS_FAIL_IF(tag != inTag, "tags do not match");

    //unpack destination part Id
    int destPartId;
    PCUObj->Unpack(destPartId);
    PCU_ALWAYS_ASSERT(rank == destPartId);

    //unpack array length
    size_t n;
    PCUObj->Unpack(n);

    //unpack int array
    int buff;
    for(size_t i=0; i < n; i++) {
      PCUObj->Unpack(buff);
      msg.push_back(buff);
    }
  }

  void recvIntsFromNeighbors(vector<int>& msg, int tag, pcu::PCU *PCUObj) {
    while(PCUObj->Listen())
      unpackInts(msg, tag, PCUObj);
  }

  int sendNetToNeighbors(partInfo& part, pcu::PCU *PCUObj) {
    PCUObj->Begin();
    MIS_ITERATE(vector<int>, part.adjPartIds, adjPartIdItr) {
      const int destRank = *adjPartIdItr;
      PCUObj->Pack(destRank, part.id); //FIXME redundant
      PCUObj->Pack(destRank, *adjPartIdItr); //FIXME redundant
      PCUObj->Pack(destRank, part.randNum);
      size_t netSz = part.net.size();
      PCUObj->Pack(destRank, netSz);
      MIS_ITERATE(vector<int>, part.net, pItr)
        PCUObj->Pack(destRank, *pItr);
    }
    return 0;
  }

  void unpackNet(vector<adjPart>& msg, pcu::PCU *PCUObj) {
    int srcPartId;
    PCUObj->Unpack(srcPartId);
    PCU_ALWAYS_ASSERT(PCUObj->Sender() == srcPartId);

    int destPartId;
    PCUObj->Unpack(destPartId);
    PCU_ALWAYS_ASSERT(PCUObj->Self() == destPartId);

    unsigned randNum;
    PCUObj->Unpack(randNum);

    adjPart ap;
    ap.partId = srcPartId;
    ap.randNum = randNum;

    size_t arrayLen;
    PCUObj->Unpack(arrayLen);
    PCU_ALWAYS_ASSERT(arrayLen > 0);

    int buff;
    for(size_t i=0; i < arrayLen; i++) {
      PCUObj->Unpack(buff);
      ap.net.push_back(buff);
    }

    msg.push_back(ap);
  }

  void recvNetsFromNeighbors(vector<adjPart>& msg, pcu::PCU *PCUObj) {
    while( PCUObj->Listen() )
      unpackNet(msg, PCUObj);
  }

  void unpackAdjPart(vector<adjPart>& msg, pcu::PCU *PCUObj) {
    const int rank = PCUObj->Self();

    //unpack destination part Id
    int destPartId;
    PCUObj->Unpack(destPartId);
    PCU_ALWAYS_ASSERT(rank == destPartId);

    //unpack adjacent part Id
    int adjPartId;
    PCUObj->Unpack(adjPartId);

    //unpack random number
    unsigned randNum;
    PCUObj->Unpack(randNum);

    adjPart ap;
    ap.partId = adjPartId;
    ap.randNum = randNum;

    //unpack net size
    size_t n;
    PCUObj->Unpack(n);

    //unpack net
    int pid;
    for(size_t i=0; i < n; i++) {
      PCUObj->Unpack(pid);
      ap.net.push_back(pid);
    }

    msg.push_back(ap);
  }

  void recvAdjNetsFromNeighbors(vector<adjPart>& msg, pcu::PCU *PCUObj) {
    while( PCUObj->Receive() )
      unpackAdjPart(msg, PCUObj);
  }

  int minRandNum(mapiuItr first, mapiuItr last) {
    unsigned min = std::numeric_limits<unsigned>::max();
    for (mapiuItr itr = first; itr != last; itr++) {
      if (itr->second < min) {
        min = itr->second;
      }
    }
    return min;
  }

  void getNetAdjPartIds(mapiuItr first, mapiuItr last, vector<int>& ids) {
    for (mapiuItr itr = first; itr != last; itr++) {
      ids.push_back(itr->first);
    }
  }
  void getNetPartIds(vector<int>& net, vector<int>&ids) {
    for (unsigned int i=0;i<net.size();i++)
      ids.push_back(net[i]);
  }

  void setNodeStateInGraph(partInfo& part) {
    if( part.adjPartIds.size() ) {
      part.isInNetGraph = true;
      part.isInMIS = false;
    } else {
      part.isInNetGraph = false;
      part.isInMIS = false;
    }
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
   * @param part (InOut) local part
   * @return 0 on success, non-zero otherwise
   */
  int constructNetGraph(partInfo& part, pcu::PCU *PCUObj) {
    stable_sort(part.adjPartIds.begin(), part.adjPartIds.end());
    setNodeStateInGraph(part);

    // send net to each neighbor
    int ierr = sendNetToNeighbors(part, PCUObj);
    if (ierr != 0) return ierr;
    PCUObj->Send();
    vector<adjPart> nbNet;
    recvNetsFromNeighbors(nbNet, PCUObj);

    // get the random num associated with adjacent nets
    part.addNetNeighbors(nbNet);

    ierr = sendAdjNetsToNeighbors(part, nbNet, PCUObj);
    if (ierr != 0) return ierr;
    PCUObj->Send();
    vector<adjPart> nbAdjNets;
    recvAdjNetsFromNeighbors(nbAdjNets, PCUObj);

    // get the random num associated with adjacent nets
    part.addNetNeighbors(nbAdjNets);

    //add new neighbors based on changes to net-graph neighbors
    part.updateNeighbors();

    return 0;
  }
}//end namespace

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
  MIS_ITERATE(mapiu, netAdjParts, netAdjPartItr) {
    findItr = find(adjPartIds.begin(),
                   adjPartIds.end(),
                   netAdjPartItr->first);
    if (findItr == adjPartIds.end())
      adjPartIds.push_back(netAdjPartItr->first);
  }
  stable_sort(adjPartIds.begin(), adjPartIds.end());
}

void mis_init(unsigned randNumSeed, int, const char*,
    const char*, const char*) {
  mersenne_twister_seed(randNumSeed);
}

void misFinalize() {}

/** If isNeighbors is true then the net of a selected part is removed from
 * the graph, o.w. the parts that have overlapping nets, as listed in
 * netAdjParts, are removed from the graph.
 *
 * Setting isNeighbors true supports computing on the partition model graph
 * as opposed to the netgraph.
 */
int mis(partInfo& part, pcu::PCU *PCUObj, bool randNumsPredefined, bool isNeighbors) {
  PCU_ALWAYS_ASSERT(PCUObj != nullptr);

  if (false == randNumsPredefined)
    setRandomNum(part, PCUObj);

  int ierr = constructNetGraph(part, PCUObj);
  if (ierr != 0) return ierr;

  vector<int> nodesRemoved;
  vector<int> nodesToRemove;
  vector<int> rmtNodesToRemove;

  int isInMis = 0;
  int tag = 0;
  int numNodesAdded;
  do {
    numNodesAdded = 0;
    const unsigned minRand =
      minRandNum(part.netAdjParts.begin(), part.netAdjParts.end());
    // add self to MIS
    if (true == part.isInNetGraph &&
        false == part.isInMIS &&
        part.randNum < minRand) {
      part.isInMIS = true;
      isInMis = 1;
      ++numNodesAdded;
      if (isNeighbors) {
        nodesToRemove.reserve(part.net.size()+1);
        getNetPartIds(part.net,nodesToRemove);
      }
      else {
        nodesToRemove.reserve(part.netAdjParts.size() + 1);
        getNetAdjPartIds(part.netAdjParts.begin(), part.netAdjParts.end(),
            nodesToRemove);
        nodesToRemove.push_back(part.id);
      }
    }

    //--------------ROUND
    tag++;
    sendIntsToNeighbors(part, nodesToRemove, tag, PCUObj);
    PCUObj->Send();
    recvIntsFromNeighbors(rmtNodesToRemove, tag, PCUObj);


    if (true == part.isInNetGraph &&
        ( true == part.isInMIS ||
          find(rmtNodesToRemove.begin(),
            rmtNodesToRemove.end(),
            part.id) != rmtNodesToRemove.end())
        ) {
      part.netAdjParts.clear();
      part.isInNetGraph = false;
      nodesRemoved.push_back(part.id);
    }
    removeNodes(part, rmtNodesToRemove);
    nodesToRemove.clear();
    rmtNodesToRemove.clear();

    //--------------ROUND

    tag++;
    sendIntsToNeighbors(part, nodesRemoved, tag, PCUObj);
    PCUObj->Send();
    recvIntsFromNeighbors(rmtNodesToRemove, tag, PCUObj);

    removeNodes(part, rmtNodesToRemove);
    nodesRemoved.clear();
    rmtNodesToRemove.clear();

    numNodesAdded = PCUObj->Add<int>(numNodesAdded);
  } while (numNodesAdded > 0);

  return isInMis;
}
