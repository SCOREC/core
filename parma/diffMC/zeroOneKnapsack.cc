// http://stackoverflow.com/questions/11036689/solving-the-integer-knapsack

#include "zeroOneKnapsack.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <math.h>

using std::cout;
using std::cin;
using std::endl;
using std::max;
using std::vector;

knapsack::knapsack(int MaxWeight, int NumItems, int* itemWeight, int* itemValue) {
    numItems = NumItems;
    maxWeight = MaxWeight;
    weight = itemWeight;
    value = itemValue;
    
    M = (int**) calloc(numItems,sizeof(int*));
    for(int i=0; i<numItems; i++) {
        M[i] = (int*) calloc(maxWeight+1,sizeof(int));
    }
}

knapsack::~knapsack() {
    for (int i = 0; i < numItems; i++) {
        free(M[i]);
    }
    free(M);
}

int knapsack::printTable() {
  printf("===== Table =====\n");
  printf("%3s | ", "");
  for(int j = 0; j <numItems; j++){
    printf("%3d", weight[j]);
  }
  printf("\n");
  printf("%3s  ","---");
  for(int j = 0; j <numItems; j++){
    printf("%3s","---");
  }
  printf("\n");
  for(int i = 1; i <= maxWeight; i++){
    printf("%3d | ", i);
    for(int j = 0; j <numItems; j++){
      printf("%3d", M[j][i]);
    }
    printf("\n");
  }
  printf("===== Table =====\n");
  return 0;
}

int knapsack::getSolution(vector<int>& solnIdx) {
    int i = maxWeight;
    int j = numItems - 1;
    while (j >= 0 && M[j][i] != 0) {
        // item 0 was taken
        if (j == 0 && weight[j] <= i && M[j][i] == value[j]) {
            solnIdx.push_back(j);
            i = i - weight[j];
            j = j - 1;
        }
        // item j!=0 was taken
        else if (weight[j] <= i && M[j][i] == M[j - 1][i - weight[j]] + value[j]) {
            solnIdx.push_back(j);
            i = i - weight[j];
            j = j - 1;
        }
        // item j was not taken
        else {
            j = j - 1;
        }
    }
    return 0;
}

int knapsack::solve() {
  for(int i = 1; i <= maxWeight; i++){
    for(int j = 0; j <numItems; j++){
      if(j > 0){
        M[j][i] = M[j-1][i];
        if (weight[j] <= i)
          M[j][i] = max(M[j][i], M[j-1][i-weight[j]]+value[j]);
      }
      else{
        M[j][i] = 0;
        if(weight[j] <= i)
          M[j][i] = max(M[j][i], value[j]);
      }
    }
  }        
    
  return M[numItems-1][maxWeight];
}  



