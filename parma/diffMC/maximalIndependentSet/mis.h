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

    typedef std::map<int, int> mapIntInt;
    typedef std::map<int, int>::iterator netAdjItr;

    typedef struct AdjPart {
        int partId;
        int randNum;
        std::vector<int> net;
    } adjPart;

    typedef struct PartInfo {
        int id;
        std::vector<int> adjPartIds;
        std::vector<int> net;

        int randNum;
        bool isInNetGraph;
        bool isInMIS;
        std::map<int, int> netAdjParts; // (partId, randNum)
        void addNetNeighbors(std::vector<adjPart>& nbNet);
        void updateNeighbors();
        void print();
    } partInfo;

    typedef std::vector<partInfo>::iterator partInfoVecItrType;

    //from http://learningcppisfun.blogspot.com/2007/02/
    //    functors-with-state-3-print-contents-of.html
    template<typename T, typename InputIterator>
    void Print(std::ostream& ostr, std::string dbgMsg, InputIterator itbegin,
        InputIterator itend, const std::string& delimiter, bool dbg=false) {
        if (dbg) {
            ostr << "DEBUG " << dbgMsg;
            std::ostream_iterator<T> out_it(ostr, delimiter.c_str());
            std::copy(itbegin, itend, out_it);
            ostr << "\n";
        }
    }
    void Print(std::ostream& ostr, char* dbgMsg,
        std::vector<misLuby::adjPart>::iterator itbegin,
        std::vector<misLuby::adjPart>::iterator itend,
        const std::string& delimiter, bool dbg=false);
    void Print(std::ostream& ostr, std::string dbgMsg,
        std::map<int, int>::iterator itbegin,
        std::map<int, int>::iterator itend,
        const std::string& delimiter, bool dbg=false);
} //end misLuby namespace

/**
 * @brief generate randNums.size() random numbers
 * @param randNums (InOut) random numbers
 * @return 0 on success, non-zero otherwise
 */
int generateRandomNumbers(std::vector<int>& randNums);

/**
 * @brief compute the maximal independent set
 * @param part (In) info on local part
 * @param randNumsPredefined (In) 0: compute random numbers, 1:uses defined random numbers
 * @return 1 if local part is in mis, 0 o.w.
 */
int mis(misLuby::partInfo& part,
    bool randNumsPredefined = false,
    bool isNeighbors = false);

void mis_init(unsigned int randNumSeed, int debugMode = 0, const char* maj = "1",
    const char* min = "0", const char* patch = "0");

void misFinalize();

#endif
