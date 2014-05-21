#include "parma_hist.h"
#include <stdio.h>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <fstream>
#include <iostream>

hist::hist() {
}

hist::~hist() {
}

/**
 * @brief add a data point
 * @param p (In) data point
 */
void hist::add(const int p) {
  d.push_back(p);
}

/**
 * @brief get the number of stored points
 * @return number of points
 */
int hist::getSz() {
  return d.size();
}

/**
 * @brief write the sorted data 
 * @param n (In) number of bins
 */
void hist::print(const int n, std::string outfile) {
  if( 0 == d.size() ) return;
  std::sort(d.begin(), d.end());

  std::ofstream ofs;
  ofs.open(outfile.c_str(), std::ofstream::out);

  std::stringstream s;
  for(size_t i=0; i<d.size(); i++) 
    s << ' ' << d[i] << ' ';
  ofs << s.str() << '\n';
  s.str("");

  const int r = d[d.size()-1] - d[0];
  const int br = ceil((double)r/(double)n);
  size_t i = 0;
  int st = d[0];
  for(int b=0; b<n; b++) {
    const int l = i;
    while( d[i] <= br*(b+1)+st && i < d.size()) 
     i++;
    const char f = b?'(':'[';
    s << f << br*b+st << ':' << br*(b+1)+st << "] " << i - l << '\n';
  }
  ofs << s.str();
}

/**
 * @brief clear added data
 */
void hist::clear() {
  d.clear();
}
