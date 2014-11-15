#include <stdio.h>
#include <set>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "mis.h"
#include "math.h"
#include <unistd.h>
#include <string>
#include "PCU.h"
#include "parma_commons.h"

using std::set;
using std::vector;
using std::find;

using std::cout;

using namespace misLuby;
using parmaCommons::debug;
using parmaCommons::status;
using parmaCommons::error;

int getPartIdFromIdx(int row, int col, const int sqrtTotNumParts) {
    if (row < 0) {
        row = sqrtTotNumParts + (row % sqrtTotNumParts);
    } else if (row >= sqrtTotNumParts) {
        row %= sqrtTotNumParts;
    }
    if (col < 0) {
        col = sqrtTotNumParts + (col % sqrtTotNumParts);
    } else if (col >= sqrtTotNumParts) {
        col %= sqrtTotNumParts;
    }
    return row * sqrtTotNumParts + col;
}

void set2dStencilNets(vector<int>& nets, const int r, const int c, 
    const int sqrtTotNumParts) {
    int id = 0;
    id = getPartIdFromIdx(r, c, sqrtTotNumParts); //self
    nets.push_back(id);
    id = getPartIdFromIdx(r, c - 1, sqrtTotNumParts); //W
    nets.push_back(id);
    id = getPartIdFromIdx(r - 1, c, sqrtTotNumParts); //N
    nets.push_back(id);
    id = getPartIdFromIdx(r, c + 1, sqrtTotNumParts); //E
    nets.push_back(id);
    id = getPartIdFromIdx(r + 1, c, sqrtTotNumParts); //S
    nets.push_back(id);
}

/**
 * @brief check that the net of part (r,c) is not in the IS
 * @param isInMIS (In) ids of parts that are in the IS
 * @param r (In) row of part being checked
 * @param c (In) column of part being checked
 * @param sqrtTotNumParts (In) floor(sqrt(totNumParts))
 * @return true if the net of part (r,c) is in the IS, false o.w.
 */
bool is2dStencilNetInIS(vector<int>& isInMIS, const int r, const int c, 
    const int sqrtTotNumParts) {
    int id = 0;
    id = getPartIdFromIdx(r, c - 1, sqrtTotNumParts); //W
    if (find(isInMIS.begin(), isInMIS.end(), id) != isInMIS.end()) {
        return true;
    }
    id = getPartIdFromIdx(r - 1, c, sqrtTotNumParts); //N
    if (find(isInMIS.begin(), isInMIS.end(), id) != isInMIS.end()) {
        return true;
    }
    id = getPartIdFromIdx(r, c + 1, sqrtTotNumParts); //E
    if (find(isInMIS.begin(), isInMIS.end(), id) != isInMIS.end()) {
        return true;
    }
    id = getPartIdFromIdx(r + 1, c, sqrtTotNumParts); //S
    if (find(isInMIS.begin(), isInMIS.end(), id) != isInMIS.end()) {
        return true;
    }

    return false;
}

void set2dStencilNeighbors(vector<int>& adjPartIds, const int r, const int c, 
    const int sqrtTotNumParts) {
    int id = 0;
    id = getPartIdFromIdx(r, c - 1, sqrtTotNumParts); //W
    adjPartIds.push_back(id);
    id = getPartIdFromIdx(r - 1, c - 1, sqrtTotNumParts); //NW
    adjPartIds.push_back(id);
    id = getPartIdFromIdx(r - 1, c, sqrtTotNumParts); //N
    adjPartIds.push_back(id);
    id = getPartIdFromIdx(r - 1, c + 1, sqrtTotNumParts); //NE
    adjPartIds.push_back(id);
    id = getPartIdFromIdx(r, c + 1, sqrtTotNumParts); //E
    adjPartIds.push_back(id);
    id = getPartIdFromIdx(r + 1, c + 1, sqrtTotNumParts); //SE
    adjPartIds.push_back(id);
    id = getPartIdFromIdx(r + 1, c, sqrtTotNumParts); //S
    adjPartIds.push_back(id);
    id = getPartIdFromIdx(r + 1, c - 1, sqrtTotNumParts); //SW
    adjPartIds.push_back(id);
}

bool isValidIndependentSet(int* isInMIS, const int sqrtTotNumParts, 
    const int totNumParts) {
    vector<int> mis;
    for (int partId = 0; partId < totNumParts; partId++) {
        if (1 == isInMIS[partId]) {
            mis.push_back(partId);
        }
    }
    for (int partId = 0; partId < totNumParts; partId++) {
        if (1 == isInMIS[partId]) {
            const int row = partId / sqrtTotNumParts;
            const int col = partId % sqrtTotNumParts;
            if (true == is2dStencilNetInIS(mis, row, col, sqrtTotNumParts)) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief determine if the net of part (r,c) is available
 * @param isInMIS (In) ids of parts that are in the IS
 * @param r (In) row of part being checked
 * @param c (In) column of part being checked
 * @param sqrtTotNumParts (In) floor(sqrt(totNumParts))
 * @return true if the net of part (r,c) is available, false o.w.
 */
bool is2dStencilNetAvailable(vector<int>& isInMIS, const int r, const int c, 
    const int sqrtTotNumParts) {
    int id = 0;
    int nb_r, nb_c;

    nb_r = r;
    nb_c = c - 1;
    id = getPartIdFromIdx(nb_r, nb_c, sqrtTotNumParts); //W
    if (find(isInMIS.begin(), isInMIS.end(), id) != isInMIS.end() ||
        true == is2dStencilNetInIS(isInMIS, nb_r, nb_c, sqrtTotNumParts)) {
      return false;
    }

    nb_r = r - 1;
    nb_c = c;
    id = getPartIdFromIdx(r - 1, c, sqrtTotNumParts); //N
    if (find(isInMIS.begin(), isInMIS.end(), id) != isInMIS.end() ||
        true == is2dStencilNetInIS(isInMIS, nb_r, nb_c, sqrtTotNumParts)) {
      return false;
    }

    nb_r = r;
    nb_c = c + 1;
    id = getPartIdFromIdx(r, c + 1, sqrtTotNumParts); //E
    if (find(isInMIS.begin(), isInMIS.end(), id) != isInMIS.end() ||
        true == is2dStencilNetInIS(isInMIS, nb_r, nb_c, sqrtTotNumParts)) {
      return false;
    }

    nb_r = r + 1;
    nb_c = c;
    id = getPartIdFromIdx(r + 1, c, sqrtTotNumParts); //S
    if (find(isInMIS.begin(), isInMIS.end(), id) != isInMIS.end() ||
        true == is2dStencilNetInIS(isInMIS, nb_r, nb_c, sqrtTotNumParts)) {
      return false;
    }

    return true;
}

bool isMaximalIndependentSet(int* isInMIS, const int sqrtTotNumParts, 
    const int totNumParts) {
  vector<int> mis;
  for (int partId = 0; partId < totNumParts; partId++) {
    if (1 == isInMIS[partId]) {
      mis.push_back(partId);
    }
  }
  for (int partId = 0; partId < totNumParts; partId++) {
    if (0 == isInMIS[partId]) {
      const int row = partId / sqrtTotNumParts;
      const int col = partId % sqrtTotNumParts;
      if (true == is2dStencilNetAvailable(mis, row, col, sqrtTotNumParts)) {
        return false;
      }
    }
  }
  return true;
}

/* 
---3x3---
 
  6 7 8   
8 _____ 6
2|0 1 2|0
5|3 4 5|3
8|6 7 8|6
2 ----- 0
  0 1 2

If any part is selected then the remaining parts are removed.
|MIS| = 1.

---4x4---

   12 13 14 15
15 ___________ 12
3 |0  1  2  3 |0
7 |4  5  6  7 |4
11|8  9  10 11|8
15|12 13 14 15|12
3  -----------
   0  1  2  3 

If part 5 is selected then first remove first adjacent neighbors in the 
netgraph:

   12 13 14 15
15 ___________ 12
3 |0     2  3 |0
7 |   5     7 | 
11|8     10 11|8
15|12 13 14 15|12
3  -----------
   0     2  3 

Now remove second adjacent neighbors of part 5:

   12    14 15
15 ___________ 12
3 |         3 | 
  |   5       | 
11|         11| 
15|12    14 15|12
3  -----------
            3 

Then suppose part 11 was selected and its first and second adjacent neighbors 
are removed:
              
   ___________   
  |           | 
  |   5       | 
11|         11| 
  |           |  
   -----------
              
|MIS| = 2

*/
int test_2dStencil(const int rank, const int totNumParts,
    const int numParts, bool randNumsPredefined = false) {
  assert(numParts==1);
  char outMsg[256];

  if( !PCU_Comm_Self() )
    status("test_2dStencil - totNumParts: %d  numParts: %d\n", 
        totNumParts, numParts);

  const int sqrtTotNumParts = floor(sqrt(totNumParts));
  if (totNumParts != sqrtTotNumParts * sqrtTotNumParts) {
    MIS_FAIL("totNumParts must be a perfect square\n");
    return 1;
  }
  if (totNumParts < 9) {
    MIS_FAIL("totNumParts must be at least 9\n");
    return 1;
  }

  //set neighbors and nets
  vector<partInfo> parts(numParts);
  partInfoVecItrType partItr;
  int partId = rank*numParts;
  for (partItr = parts.begin(); partItr != parts.end(); partItr++) {
    const int row = partId / sqrtTotNumParts;
    const int col = partId % sqrtTotNumParts;
    partItr->id = partId;
    set2dStencilNeighbors(partItr->adjPartIds, row, col, sqrtTotNumParts);
    set2dStencilNets(partItr->net, row, col, sqrtTotNumParts);
    partId++;
  }

  //sanity check
  for (partItr = parts.begin(); partItr != parts.end(); partItr++) {
    vector<int>::iterator adjPIDItr;
    for (adjPIDItr = partItr->adjPartIds.begin();
        adjPIDItr != partItr->adjPartIds.end();
        adjPIDItr++) {
      if (*adjPIDItr < 0 || *adjPIDItr >= totNumParts) {
        printf("ERROR [%d] (%d) adjPartId=%d\n", 
            rank, partItr->id, *adjPIDItr);
        return 1;
      }
    }
    vector<int>::iterator netItr;
    for (netItr = partItr->net.begin();
        netItr != partItr->net.end();
        netItr++) {
      if (*netItr < 0 || *netItr >= totNumParts) {
        printf("ERROR [%d] (%d) net=%d\n", rank, partItr->id, *netItr);
        return 1;
      }
    }
  }

  vector<int> maximalIndSet;
  const double t1 = MPI_Wtime();
  int ierr = mis(parts[0], maximalIndSet, randNumsPredefined);
  double elapsedTime = MPI_Wtime() - t1;
  PCU_Max_Doubles(&elapsedTime, 1);

  if( !PCU_Comm_Self() )
    status("elapsed time (seconds) = %f \n", elapsedTime);

  vector<int> isInMIS(numParts, 0);
  vector<int>::iterator misItr;
  for (misItr = maximalIndSet.begin();
      misItr != maximalIndSet.end();
      misItr++) {
    const int partIdx = *misItr % numParts;
    isInMIS[partIdx] = 1;
  }

  int* globalIsInMIS = NULL;
  if (rank == 0) {
    globalIsInMIS = new int[totNumParts];
  }
  MPI_Gather(&(isInMIS[0]), numParts, MPI_INT, globalIsInMIS, numParts, 
      MPI_INT, 0, MPI_COMM_WORLD);
  int sizeIS = 0;

  if (rank != 0) 
    return 0;

  char* statusMsg = new char[totNumParts * 4];
  int pos = sprintf(statusMsg, "MIS: ");
  for (int i = 0; i < totNumParts; ++i) {
    if (globalIsInMIS[i]) {
      pos += sprintf(&(statusMsg[pos]), " %d ", i);
      ++sizeIS;
    }
  }
  status("%s\n", statusMsg);
  delete [] statusMsg;

  status("|independent set| = %d\n", sizeIS);

  if (!isValidIndependentSet(globalIsInMIS, sqrtTotNumParts, totNumParts))
    MIS_FAIL("Not a valid independent set!\n");
  if (!isMaximalIndependentSet(globalIsInMIS, sqrtTotNumParts, totNumParts))
    MIS_FAIL("Not a maximal independent set!\n");

  delete [] globalIsInMIS;

  int error = 1;
  if ( totNumParts == 9 && sizeIS == 1 )
    return !error;
  else if ( totNumParts == 16 && sizeIS == 2 )
    return !error;
  else if ( totNumParts > 16 && sizeIS > 2 )
    return !error;
  else 
    return error;
}

/**
 * @brief part 0, on rank 0, is the center of the star, all other ranks are satellites
 * @remark net(sattelite i) = {0, i}
 *         neighbors(sattelite i) = {0, i+1, i-1}
 *         net(0) = {all satellites}
 *         neighbors(0) = {all satellites}
 *
 *         the size of the MIS = 1
 * 
 * @param rank (In) MPI process rank
 * @param totNumParts (In) total number of parts
 * @param numParts (In) number of local parts
 * @param randNumsPredefined (In) T:rand nums predefined, F:rand nums not predefined
 * @return 0 on success, non-zero otherwise
 */
int test_StarA(const int rank, const int totNumParts,
        const int numParts, bool randNumsPredefined = false) {
    char dbgMsg[256];
    assert(numParts == 1);

    if( !PCU_Comm_Self() )
      status("test_StarA - totNumParts: %d  numParts: %d\n", totNumParts, 
          numParts);

    vector<partInfo> parts;

    for (int partIdx = 0; partIdx < numParts; partIdx++) {
        partInfo p;   
        if (0 == rank && 0 == partIdx) {
            p.id = 0;
            p.net.push_back(0);
            for (int partId = 1; partId < totNumParts; ++partId) {
                p.adjPartIds.push_back(partId);
                p.net.push_back(partId);
            }
        } else {
            p.id = rank * numParts + partIdx;

            p.net.push_back(p.id);
            p.adjPartIds.push_back(0);

            if (1 == p.id) {
                p.adjPartIds.push_back(totNumParts - 1);
            } else {
                p.adjPartIds.push_back(p.id - 1);
            }

            if (totNumParts - 1 == p.id) {
                p.adjPartIds.push_back(1);
            } else {
                p.adjPartIds.push_back(p.id + 1);
            }
        }
        parts.push_back(p);
    }
    
    assert(parts.size() == numParts);

    partInfoVecItrType partItr;
    for (partItr = parts.begin(); partItr != parts.end(); partItr++) {
        vector<int>::iterator adjPIDItr;
        for (adjPIDItr = partItr->adjPartIds.begin();
                adjPIDItr != partItr->adjPartIds.end();
                adjPIDItr++) {
            if (*adjPIDItr < 0 || *adjPIDItr >= totNumParts) {
                printf("ERROR [%d] adjPartId=%d\n", rank, *adjPIDItr);
                return 1;
            }
        }
        sprintf(dbgMsg, " [%d] (%d) %s - adjPartIds=", 
            rank, partItr->id, __FUNCTION__);
        Print <int, vector<int>::iterator > (cout, dbgMsg, 
            partItr->adjPartIds.begin(), partItr->adjPartIds.end(), 
            std::string(", "));

        sprintf(dbgMsg, " [%d] (%d) %s - net=", 
            rank, partItr->id, __FUNCTION__);
        Print <int, vector<int>::iterator > (cout, dbgMsg, 
            partItr->net.begin(), partItr->net.end(), std::string(", "));
    }

    vector<int> maximalIndSet;
    int ierr = mis(parts[0], maximalIndSet, randNumsPredefined);

    vector<int> isInMIS(numParts, 0);
    vector<int>::iterator misItr;
    for (misItr = maximalIndSet.begin();
            misItr != maximalIndSet.end();
            misItr++) {
        const int partIdx = *misItr % numParts;
        isInMIS[partIdx] = 1;
    }

    int* globalIsInMIS = NULL;
    if (rank == 0) {
        globalIsInMIS = new int[totNumParts];
    }
    MPI_Gather(&(isInMIS[0]), numParts, MPI_INT, globalIsInMIS, numParts, 
        MPI_INT, 0, MPI_COMM_WORLD);
    int sizeIS = 0;
    if (rank == 0) {
        char* statusMsg = new char[totNumParts * 4];
        int pos = sprintf(statusMsg, "MIS: ");
        for (int i = 0; i < totNumParts; ++i) {
            if (globalIsInMIS[i]) {
                pos += sprintf(&(statusMsg[pos]), " %d ", i);
                ++sizeIS;
            }
        }
        status("%s\n", statusMsg);
        delete [] globalIsInMIS;

        status("|independent set| = %d\n", sizeIS);
        delete [] statusMsg;


        if (sizeIS == totNumParts - 1 || sizeIS == 1) {
            return 0;
        } else {
            return 1;
        }
    } else {
        return 0;
    }
}

/*

part and net graph
    0 1 2 3 
  0   x x x
  1 x     x
  2 x     x
  3 x x x

part      0  1  2  3
randNums  5  1  13 6

  0 1 2 3
  -------
    
 */
int test_4partsA(const int rank, const int totNumParts, const int numParts, 
    bool predefinedRandNums) {
    assert(numParts == 1);
    if (4 != totNumParts) {
        MIS_FAIL("totNumParts must be 4\n");
        return 1;
    }

    char dbgMsg[256];
    if( !PCU_Comm_Self() )
      status("test_4partsA - totNumParts: %d  numParts: %d\n", 
          totNumParts, numParts);

    partInfo part;
    part.id = rank;

    if (rank == 0) {
        part.randNum = 5;
        part.net.push_back(0);
        part.net.push_back(1);
        part.net.push_back(2);
        part.adjPartIds.push_back(1);
        part.adjPartIds.push_back(2);
        part.adjPartIds.push_back(3);
    } else if (rank == 1) {
        part.randNum = 1;
        part.net.push_back(0);
        part.net.push_back(1);
        part.adjPartIds.push_back(0);
        part.adjPartIds.push_back(3);
    } else if (rank == 2) {
        part.randNum = 13;
        part.net.push_back(2);
        part.net.push_back(3);
        part.adjPartIds.push_back(0);
        part.adjPartIds.push_back(3);
    } else if (rank == 3) {
        part.randNum = 6;
        part.net.push_back(2);
        part.net.push_back(3);
        part.adjPartIds.push_back(0);
        part.adjPartIds.push_back(1);
        part.adjPartIds.push_back(2);
    }

    part.print();

    vector<partInfo> parts;
    parts.push_back(part);
    vector<int> maximalIndSet;
    int ierr = mis(parts[0], maximalIndSet, predefinedRandNums);

    int misLocalSize = maximalIndSet.size();
    int* globalIsInMIS = NULL;
    if (rank == 0) {
        globalIsInMIS = new int[totNumParts];
    }
    MPI_Gather(&misLocalSize, 1, MPI_INT, globalIsInMIS, 1, MPI_INT, 
        0, MPI_COMM_WORLD);
    int sizeIS = 0;
    if (rank == 0) {
        printf("MIS: ");
        for (int i = 0; i < totNumParts; ++i) {
            if (globalIsInMIS[i]) {
                printf(" %d ", i);
                ++sizeIS;
            }
        }
        printf("\n");
        delete [] globalIsInMIS;

        if (sizeIS < totNumParts / 3) {
            return 1;
        } else {
            printf("\n\nDEBUG |independent set| = %d\n", sizeIS);
            return 0;
        }
    } else {
        return 0;
    }
}

int test_4partsB(const int rank, const int totNumParts, const int numParts, 
    bool predefinedRandNums) {
    if (4 != totNumParts) {
        MIS_FAIL("totNumParts must be 4\n");
        return 1;
    }
    assert(numParts==1);

    char dbgMsg[256];
    if( !PCU_Comm_Self() )
      status("test_4partsAndNums - totNumParts: %d  numParts: %d\n", 
          totNumParts, numParts);

    partInfo part;
    part.id = rank;

    if (rank == 0) {
        part.randNum = 1;
        part.net.push_back(0);
        part.adjPartIds.push_back(1);
        part.adjPartIds.push_back(2);
        part.adjPartIds.push_back(3);
    } else if (rank == 1) {
        part.randNum = 4;
        part.net.push_back(0);
        part.net.push_back(1);
        part.adjPartIds.push_back(0);
        part.adjPartIds.push_back(3);
    } else if (rank == 2) {
        part.randNum = 13;
        part.net.push_back(2);
        part.net.push_back(0);
        part.adjPartIds.push_back(0);
        part.adjPartIds.push_back(3);
    } else if (rank == 3) {
        part.randNum = 6;
        part.net.push_back(0);
        part.net.push_back(3);
        part.adjPartIds.push_back(0);
        part.adjPartIds.push_back(1);
        part.adjPartIds.push_back(2);
    }

    vector<partInfo> parts;
    parts.push_back(part);
    vector<int> maximalIndSet;
    int ierr = mis(parts[0], maximalIndSet, predefinedRandNums);

    int misLocalSize = maximalIndSet.size();
    int* globalIsInMIS = NULL;
    if (rank == 0) {
        globalIsInMIS = new int[totNumParts];
    }
    MPI_Gather(&misLocalSize, 1, MPI_INT, globalIsInMIS, 1, MPI_INT, 
        0, MPI_COMM_WORLD);
    int sizeIS = 0;
    if (rank == 0) {
        printf("MIS: ");
        for (int i = 0; i < totNumParts; ++i) {
            if (globalIsInMIS[i]) {
                printf(" %d ", i);
                ++sizeIS;
            }
        }
        printf("\n");
        delete [] globalIsInMIS;

        if (sizeIS < totNumParts / 3) {
            return 1;
        } else {
            printf("\n\nDEBUG |independent set| = %d\n", sizeIS);
            return 0;
        }
    } else {
        return 0;
    }
}

int test_4partsC(const int rank, const int totNumParts, const int numParts, 
    bool predefinedRandNums) {
    assert(numParts==1);
    if (4 != totNumParts) {
        MIS_FAIL("totNumParts must be 4\n");
        return 1;
    }

    char dbgMsg[256];
    if( !PCU_Comm_Self() )
      status("test_4partsA - totNumParts: %d  numParts: %d\n", 
          totNumParts, numParts);

    partInfo part;
    part.id = rank;

    if (rank == 0) {
        part.randNum = 1;
        part.net.push_back(0);
        part.adjPartIds.push_back(1);
        part.adjPartIds.push_back(2);
    } else if (rank == 1) {
        part.randNum = 2;
        part.net.push_back(1);
        part.adjPartIds.push_back(0);
        part.adjPartIds.push_back(3);
    } else if (rank == 2) {
        part.randNum = 13;
        part.net.push_back(2);
        part.net.push_back(0);
        part.adjPartIds.push_back(0);
    } else if (rank == 3) {
        part.randNum = 6;
        part.net.push_back(1);
        part.net.push_back(3);
        part.adjPartIds.push_back(1);
    }

    vector<partInfo> parts;
    parts.push_back(part);
    vector<int> maximalIndSet;
    int ierr = mis(parts[0], maximalIndSet, predefinedRandNums);

    int misLocalSize = maximalIndSet.size();
    int* globalIsInMIS = NULL;
    if (rank == 0) {
        globalIsInMIS = new int[totNumParts];
    }
    MPI_Gather(&misLocalSize, 1, MPI_INT, globalIsInMIS, 1, MPI_INT, 
        0, MPI_COMM_WORLD);
    int sizeIS = 0;
    if (rank == 0) {
        printf("MIS: ");
        for (int i = 0; i < totNumParts; ++i) {
            if (globalIsInMIS[i]) {
                printf(" %d ", i);
                ++sizeIS;
            }
        }
        printf("\n");
        delete [] globalIsInMIS;

        if (sizeIS < totNumParts / 3) {
            return 1;
        } else {
            printf("\n\nDEBUG |independent set| = %d\n", sizeIS);
            return 0;
        }
    } else {
        return 0;
    }
}

void printUsage(char* exe) {
    fprintf(stderr, "Usage: %s -t test [-n number of parts] [-s random number seed] [-d] \n", exe);
    fprintf(stderr, "-d enables debug mode\n"
            "-s specify an unsigned integer to seed the random number generator\n"
            "-n specify the number of parts per process, default is 1\n"
            "-r disable pre-defined random numbers for test 0, 1 and 2\n"
            "test=0 - predefined graph with 4 parts - with pre-defined random numbers MIS = {1,3}, w/random numbers MIS = {0} | {1,3} | {1,2} \n"
            "test=1 - predefined graph with 4 parts - first self net test - with pre-defined random numbers MIS = {0}, w/random numbers MIS != {0}\n"
            "test=2 - predefined graph with 4 parts - second self net test - with pre-defined random numbers MIS = {0,1}, w/random numbers MIS = {2,3}\n"
            "test=3 - part 0, on rank 0, is the center of the star, all other ranks are satellites, all satellites have a net that includes the center of the star, the MIS must be equal to one\n"
            "test=4 - each part has neighbors W,NW,N,NE,E,SE,S,SW and net W,N,E,S\n");
}

int broadcastInt(int val) {
    MPI_Bcast(&val, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return val;
}

bool getBool(int& isDebug) {
    MPI_Bcast(&isDebug, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (isDebug == 0) {
        return false;
    } else {
        return true;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    PCU_Comm_Init();
    PCU_Protect();
    int rank;
    PCU_Comm_Rank(&rank);
    int commSize;
    PCU_Comm_Size(&commSize);

    int opt;
    int debugMode = 0;
    int randNumSeed = time(NULL);
    int testNum = -1;
    int numParts = 1;
    int iPredefinedRandNums = 1;

    if (0 == rank) {
        while ((opt = getopt(argc, argv, "ds:t:n:r")) != -1) {
            switch (opt) {
                case 'd':
                    debugMode = 1;
                    break;
                case 's':
                    randNumSeed = atoi(optarg);
                    break;
                case 't':
                    testNum = atoi(optarg);
                    break;
                 case 'r':
                    iPredefinedRandNums = 0;
                    break;
                case 'n':
                    numParts = atoi(optarg);
                    break;
                default:
                    printUsage(argv[0]);
                    exit(EXIT_FAILURE);
            }
        }
    }


    numParts = broadcastInt(numParts);
    testNum = broadcastInt(testNum);
    randNumSeed = broadcastInt(randNumSeed) + PCU_Comm_Self();
    bool predefinedRandNums = getBool(iPredefinedRandNums);

    mis_init(randNumSeed, getBool(debugMode));

    int ierr;
    switch (testNum) {
      default:
      case -1:
        if (0 == rank) {
          printf("Test number not recognized\n");
          printUsage(argv[0]);
        }
        MPI_Finalize();
        return 0;
        break;
      case 0:
        ierr = test_4partsA(rank, commSize*numParts, numParts, 
            predefinedRandNums);
        break;
      case 1:
        ierr = test_4partsB(rank, commSize*numParts, numParts, 
            predefinedRandNums);
        break;
      case 2:
        ierr = test_4partsC(rank, commSize*numParts, numParts, 
            predefinedRandNums);
        break;
      case 3:
        ierr = test_StarA(rank, commSize*numParts, numParts);
        break;
      case 4:
        ierr = test_2dStencil(rank, commSize*numParts, numParts);
        break;
    }

    if (0 == rank) {
        if (1 == ierr) {
            printf("failed\n");
        } else {
            printf("passed.\n");
        }
    }
    misFinalize();
    PCU_Comm_Free();
    MPI_Finalize();
    return 0;
}
