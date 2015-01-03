#include "zeroOneKnapsack.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

int test1() {
  size_t maxw = 4;
  size_t n = 5;
  size_t w[4] = {1,2,3,4};
  size_t v[4] = {0,1,4,2};
  Knapsack k = makeKnapsack(maxw, n, w, v);
  size_t val = solve(k);
  size_t size;
  size_t* soln;
  printf("val %lu\n", val);
  assert(val == 4);
  soln = getSolution(k, &size);
  printf("size %lu soln[0] %lu\n", size, soln[0]);
  assert(size == 1);
  assert(soln[0] == 2);
  destroyKnapsack(k);
  free(soln);
  return 1;
}

int test2() {
  /* http://cse.unl.edu/~goddard/Courses/CSCE310J/Lectures/Lecture8-DynamicProgramming.pdf */
  size_t maxw = 5;
  size_t n = 4;
  size_t w[4] = {2,3,4,5};
  size_t v[4] = {3,4,5,6};
  Knapsack k = makeKnapsack(maxw, n, w, v);
  size_t val = solve(k);
  size_t size;
  size_t* soln;
  printTable(k);
  printf("val %lu\n", val);
  assert(val == 7);
  soln = getSolution(k, &size);
  printf("size %lu soln[0] %lu soln[1] %lu\n", size, soln[0], soln[1]);
  assert(size == 2);
  assert(soln[0] == 1 && soln[1] == 0);
  destroyKnapsack(k);
  free(soln);
  return 1;
}

int main() {
  if( !test1() ) return 1;
  if( !test2() ) return 1;
  return 0;
}
