#include "zeroOneKnapsack.h"
#include <stdio.h>
#include <stdlib.h>

struct zeroOneKnapsack {
    size_t **M;
    size_t *weight;
    size_t *value;
    size_t maxWeight;
    size_t numItems;
};
typedef struct zeroOneKnapsack* zoks;

size_t max(const size_t a, const size_t b);
size_t solnSz(zoks k);

size_t max(const size_t a, const size_t b) {
  if( a > b )
    return a;
  else 
    return b;
}

Knapsack makeKnapsack(size_t MaxWeight, size_t NumItems, 
    size_t* itemWeight, size_t* itemValue) {
    zoks k = (zoks) calloc(1, sizeof(struct zeroOneKnapsack));
    k->numItems = NumItems;
    k->maxWeight = MaxWeight;
    k->weight = itemWeight;
    k->value = itemValue;
    
    k->M = (size_t**) calloc(k->numItems,sizeof(size_t*));
    for(size_t i=0; i<k->numItems; i++) {
        k->M[i] = (size_t*) calloc(k->maxWeight+1,sizeof(size_t));
    }
    return (Knapsack) k;
}

void destroyKnapsack(Knapsack knapsack) {
  zoks k = (zoks) knapsack;
  for (size_t i = 0; i < k->numItems; i++) {
    free(k->M[i]);
  }
  free(k->M);
}

void printTable(Knapsack knapsack) {
  zoks k = (zoks) knapsack;
  printf("===== Table =====\n");
  printf("%3s | ", "");
  for(size_t j = 0; j <k->numItems; j++){
    printf("%3lu", k->weight[j]);
  }
  printf("\n");
  printf("%3s  ","---");
  for(size_t j = 0; j <k->numItems; j++){
    printf("%3s","---");
  }
  printf("\n");
  for(size_t i = 1; i <= k->maxWeight; i++){
    printf("%3lu | ", i);
    for(size_t j = 0; j <k->numItems; j++){
      printf("%3lu", k->M[j][i]);
    }
    printf("\n");
  }
  printf("===== Table =====\n");
}

size_t solnValue(zoks k) {
  return k->M[k->numItems-1][k->maxWeight];
}

size_t* getSolution(Knapsack knapsack, size_t* sz) {
    zoks k = (zoks) knapsack;
    size_t* soln = (size_t*) calloc(k->numItems, sizeof(size_t));
    size_t solnIdx = 0;
    size_t i = k->maxWeight;
    size_t j = k->numItems - 1;
    while ( k->M[j][i] != 0) {
        // item 0 was taken
        if (j == 0 && k->weight[j] <= i && k->M[j][i] == k->value[j]) {
            soln[solnIdx++] = j;
            i = i - k->weight[j];
        }
        // item j!=0 was taken
        else if (k->weight[j] <= i && 
            k->M[j][i] == k->M[j - 1][i - k->weight[j]] + k->value[j]) {
            soln[solnIdx++] = j;
            i = i - k->weight[j];
        }
        // item j was not taken
        if( (int)j-1 < 0 ) 
          break;
        else 
          j--;
    }
    *sz = solnIdx;
    return soln;
}

size_t solve(Knapsack knapsack) {
  zoks k = (zoks) knapsack;
  for(size_t i = 1; i <= k->maxWeight; i++){
    for(size_t j = 0; j <k->numItems; j++){
      if(j > 0){
        k->M[j][i] = k->M[j-1][i];
        if (k->weight[j] <= i)
          k->M[j][i] = max(k->M[j][i], k->M[j-1][i-k->weight[j]]+k->value[j]);
      }
      else{
        k->M[j][i] = 0;
        if(k->weight[j] <= i)
          k->M[j][i] = max(k->M[j][i], k->value[j]);
      }
    }
  }        
    
  return solnValue(k);
}  



