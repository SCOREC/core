#include "parma_commons.h"
#include <lionPrint.h>
#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <cstdarg>

///////////////////////////////////////////////////////////////////////////////////////////////
// http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/ 
inline bool AlmostEqualRelativeAndAbs(double A, double B,
            double maxDiff, double maxRelDiff)
{
    // Check if the numbers are really close -- needed
    // when comparing numbers near zero.
    double diff = fabs(A - B);
    if (diff <= maxDiff)
        return true;
 
    A = fabs(A);
    B = fabs(B);
    double largest = (B > A) ? B : A;
 
    if (diff <= largest * maxRelDiff)
        return true;
    return false;
}
///////////////////////////////////////////////////////////////////////////////////////////////   

int parmaCommons::isEql(double a, double b) {
   if( ! AlmostEqualRelativeAndAbs(a, b, 1e-8, 1e-5) )
      return 0; 
   return 1;
}

int parmaCommons::isLess(double a, double b) { 
   if( ! isEql(a, b) && a < b ) 
      return 1; 
   return 0;
}

int parmaCommons::isMore(double a, double b) { 
   if( ! isEql(a, b) && a > b ) 
      return 1; 
   return 0;
}

void parmaCommons::printElapsedTime(const char* fn, double elapsed, pcu::PCU *PCUObj) {
   elapsed = PCUObj->Max<double>(elapsed);
   if( !PCUObj->Self() )
      status("%s elapsed time %lf seconds\n", fn, elapsed);
}

void parmaCommons::debug(bool isActive, const char* fmt,...) {
  if(!isActive) return;
  lion_oprint(2,"PARMA_DEBUG ");
  lion_oprint(2,fmt);
}

void parmaCommons::status(const char* fmt,...) {
  lion_oprint(1, "PARMA_STATUS ");
  va_list ap;
  va_start(ap,fmt);
  lion_voprint(1, fmt, ap);
  va_end(ap);
}

void parmaCommons::error(const char* fmt,...) {
  lion_oprint(1,"PARMA_ERROR ");
  va_list ap;
  va_start(ap,fmt);
  lion_voprint(1, fmt, ap);
  va_end(ap);
}

