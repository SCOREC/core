#ifndef PARMA_COMMONS_H_
#define PARMA_COMMONS_H_

#include "apfArray.h"

namespace parmaCommons {

int isEql(double a, double b);
template<class type1, class type2> int isEqlArr(type1 a, type2 b) {
   for(size_t i=0; i<4; i++) 
      if( ! isEql(a[i], b[i]) )
	 return 0; 
   return 1;
}

int isLess(double a, double b);
int isMore(double a, double b);

void printElapsedTime(const char* fn, double elapsed);

void debug(bool isActive, const char* fmt,...);
void status(const char* fmt,...);
void error(const char* fmt,...);

}//end parmaCommons

#endif

