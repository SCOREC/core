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
  int generateRandomNumber() {
    return rand();
  }

  void setRandomNum(partInfo& part) {
    part.randNum = generateRandomNumber();
    // don't select self nets until all other nets are selected
    if ( 1 == part.net.size() )
      part.randNum = std::numeric_limits<int>::max();
  }

  void removeNodes(partInfo& p, vector<int>& nodes) {
    if (true == p.isInNetGraph)
      MIS_ITERATE(vector<int>, nodes, nodeItr)
        p.netAdjParts.erase(*nodeItr);
  }

  int sendAdjNetsToNeighbors(partInfo& part, vector<adjPart>& nbNet) {
    PCU_Comm_Begin();

    MIS_ITERATE(vector<int>, part.adjPartIds, adjPartIdItr) {
      const int destRank = *adjPartIdItr;
      MIS_ITERATE(vector<adjPart>, nbNet, nbItr) {
        if (nbItr->partId != *adjPartIdItr) {
          //pack destination part Id
          PCU_COMM_PACK(destRank, *adjPartIdItr);

          //adjacent part Id
          PCU_COMM_PACK(destRank, nbItr->partId);

          //adjacent part's random number
          PCU_COMM_PACK(destRank, nbItr->randNum);

          //net size
          size_t n = nbItr->net.size();
          PCU_COMM_PACK(destRank, n);

          //net
          for (vector<int>::iterator pItr = nbItr->net.begin();
              pItr != nbItr->net.end();
              pItr++) {
            PCU_COMM_PACK( destRank, *pItr);
          }
        }
      }
    }
    return 0;
  }

  void sendIntsToNeighbors(partInfo& part, vector<int>& msg, int tag) {
    PCU_Comm_Begin();
    MIS_ITERATE(vector<int>, part.adjPartIds, adjPartIdItr) {
      const int destRank = *adjPartIdItr;

      //pack msg tag
      PCU_COMM_PACK(destRank, tag);

      //pack destination part Id
      PCU_COMM_PACK(destRank, *adjPartIdItr);

      //pack array length
      size_t n = msg.size();
      PCU_COMM_PACK(destRank, n);

      //pack int array
      for ( vector<int>::iterator pItr = msg.begin();
          pItr != msg.end();
          pItr++ ) {
        PCU_COMM_PACK(destRank, *pItr);
      }
    }
  }

  void unpackInts(vector<int>& msg, int tag) {
    const int rank = PCU_Comm_Self();

    //unpack msg tag
    int inTag;
    PCU_COMM_UNPACK(inTag);
    MIS_FAIL_IF(tag != inTag, "tags do not match");

    //unpack destination part Id
    int destPartId;
    PCU_COMM_UNPACK(destPartId);
    assert(rank == destPartId);

    //unpack array length
    size_t n;
    PCU_COMM_UNPACK(n);

    //unpack int array
    int buff;
    for(size_t i=0; i < n; i++) {
      PCU_COMM_UNPACK(buff);
      msg.push_back(buff);
    }
  }

  void recvIntsFromNeighbors(vector<int>& msg, int tag) {
    while(PCU_Comm_Listen())
      unpackInts(msg, tag);
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

  void unpackNet(vector<adjPart>& msg) {
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

  void recvNetsFromNeighbors(vector<adjPart>& msg) {
    while( PCU_Comm_Listen() )
      unpackNet(msg);
  }

  void unpackAdjPart(vector<adjPart>& msg) {
    const int rank = PCU_Comm_Self();

    //unpack destination part Id
    int destPartId;
    PCU_COMM_UNPACK(destPartId);
    assert(rank == destPartId);

    //unpack adjacent part Id
    int adjPartId;
    PCU_COMM_UNPACK(adjPartId);

    //unpack random number
    int randNum;
    PCU_COMM_UNPACK(randNum);

    adjPart ap;
    ap.partId = adjPartId;
    ap.randNum = randNum;

    //unpack net size
    size_t n;
    PCU_COMM_UNPACK(n);

    //unpack net
    int pid;
    for(size_t i=0; i < n; i++) {
      PCU_COMM_UNPACK(pid);
      ap.net.push_back(pid);
    }

    msg.push_back(ap);
  }

  void recvAdjNetsFromNeighbors(vector<adjPart>& msg) {
    while( PCU_Comm_Receive() )
      unpackAdjPart(msg);
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
  int constructNetGraph(partInfo& part) {
    stable_sort(part.adjPartIds.begin(), part.adjPartIds.end());
    setNodeStateInGraph(part);

    // send net to each neighbor
    int ierr = sendNetToNeighbors(part);
    if (ierr != 0) return ierr;
    PCU_Comm_Send();
    vector<adjPart> nbNet;
    recvNetsFromNeighbors(nbNet);

    // get the random num associated with adjacent nets
    part.addNetNeighbors(nbNet);

    ierr = sendAdjNetsToNeighbors(part, nbNet);
    if (ierr != 0) return ierr;
    PCU_Comm_Send();
    vector<adjPart> nbAdjNets;
    recvAdjNetsFromNeighbors(nbAdjNets);

    // get the random num associated with adjacent nets
    part.addNetNeighbors(nbAdjNets);

    //add new neighbors based on changes to net-graph neighbors
    part.updateNeighbors();

    return 0;
  }
}//end namespace

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
  out << "randNum " << randNum << " isInNetGraph " << isInNetGraph
      << " isInMIS " << isInMIS << "\n";
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

void mis_init(unsigned int randNumSeed, int, const char*,
    const char*, const char*) {
  srand(randNumSeed);
}

void misFinalize() {}

/** If isNeighbors is true then the net of a selected part is removed from
 * the graph, o.w. the parts that have overlapping nets, as listed in
 * netAdjParts, are removed from the graph.
 *
 * Setting isNeighbors true supports computing on the partition model graph
 * as opposed to the netgraph.
 */
int mis(partInfo& part, bool randNumsPredefined,bool isNeighbors) {
  assert(PCU_Comm_Initialized());

  if (false == randNumsPredefined)
    setRandomNum(part);

  int ierr = constructNetGraph(part);
  if (ierr != 0) return ierr;

  vector<int> nodesRemoved;
  vector<int> nodesToRemove;
  vector<int> rmtNodesToRemove;

  int isInMis = 0;
  int loopCount = 0;
  int tag = 0;
  int numNodesAdded;
  do {
    numNodesAdded = 0;
    const int minRand =
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
    sendIntsToNeighbors(part, nodesToRemove, tag);
    PCU_Comm_Send();
    recvIntsFromNeighbors(rmtNodesToRemove, tag);


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
    sendIntsToNeighbors(part, nodesRemoved, tag);
    PCU_Comm_Send();
    recvIntsFromNeighbors(rmtNodesToRemove, tag);

    removeNodes(part, rmtNodesToRemove);
    nodesRemoved.clear();
    rmtNodesToRemove.clear();

    PCU_Add_Ints(&numNodesAdded, 1);
    loopCount++;
  } while (numNodesAdded > 0);

  return isInMis;
}
