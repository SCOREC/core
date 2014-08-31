#include "parma_priority.h"
#include "parma_commons.h"
#include "PCU.h"
#include <stdio.h>
#include <sstream>

using parmaCommons::status;

/// sorted entity priority entDim[i] has priority[i]
/// priority[i] > priority[i+1]
/// if priority[i] < 0 then balance of entDim[i] is not improved
void priorityList::print() {
   std::stringstream msg;
   msg << '[' << PCU_Comm_Self() << "] priorityList <entDim,priority>: ";
   for(int i=0; i<4; i++) 
      msg << '(' << entDim[i] << ',' << priority[i] << ") ";
   msg << '\n';
   std::string s = msg.str();
   status(s.c_str());
}

/**
 * @brief in-place sort of priority list
 *
 * @param priority (In) priority of entity type [vtx, edge, face, rgn]
 */
void priorityList::sort(int (*inPrio)[4]) {
   int prio[4] = {(*inPrio)[0], (*inPrio)[1], (*inPrio)[2], (*inPrio)[3]};
   
   for(int i=0; i<4; i++){
      int index = -1;
      int maxpriority = -1;
      //find the index with the maximum priority
      for(int j=0; j<4; j++){
	 if(prio[j] > maxpriority) { 
	    maxpriority = prio[j];
	    index = j;
	 } 
	 if(prio[j] == maxpriority && j > index) 
	    index = j;
      }
      assert(index != -1);
      //priority[index] has maximum priority equal to maxpriority
      priority[i] = prio[index];
      entDim[i] = index;
      prio[index] = -2;
   }
}
