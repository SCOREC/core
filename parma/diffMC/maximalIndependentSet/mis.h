#ifndef MIS_H_
#define MIS_H_

#include <map>
#include <vector>
#include <iostream>
#include <iterator>
#include <limits.h>

#include "mpi.h"

#include <stddef.h>
#include "PCU.h"

#define MIS_ITERATE(t,w,i) \
for (t::iterator (i) = (w).begin(); \
     (i) != (w).end(); ++(i))

#define MIS_FAIL(message)\
{fprintf(stderr,"MIS ERROR: %s: " message "\n", __func__);\
abort();}
#define MIS_FAIL_IF(condition,message)\
if (condition)\
MIS_FAIL(message)

namespace misLuby {
    typedef std::map<int, unsigned> mapiu;
    typedef mapiu::iterator mapiuItr;
    typedef struct AdjPart {
        int partId;
        unsigned randNum;
        std::vector<int> net;
    } adjPart;

    typedef struct PartInfo {
        int id;
        std::vector<int> adjPartIds;
        std::vector<int> net;

        unsigned randNum;
        bool isInNetGraph;
        bool isInMIS;
        mapiu netAdjParts; // (partId, randNum)
        void addNetNeighbors(std::vector<adjPart>& nbNet);
        void updateNeighbors();
    } partInfo;
} //end misLuby namespace

/**
 * @brief compute the maximal independent set
 * @param part (In) info on local part
 * @param randNumsPredefined (In) 0: compute random numbers, 1:uses defined random numbers
 * @return 1 if local part is in mis, 0 o.w.
 */
int mis(misLuby::partInfo& part,
    bool randNumsPredefined = false,
    bool isNeighbors = false);

void mis_init(unsigned randNumSeed, int debugMode = 0, const char* maj = "1",
    const char* min = "0", const char* patch = "0");

void misFinalize();

#endif
