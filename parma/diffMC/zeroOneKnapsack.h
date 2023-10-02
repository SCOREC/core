#ifndef KNAPSACK_H_
#define KNAPSACK_H_

#ifdef __cplusplus
#include <cstddef>
extern "C" {
#else
#include <stddef.h>
#endif

typedef void* Knapsack;

/**
 * @brief solve the zero-one knapsack problem
 * @remark code based on 
 *   http://stackoverflow.com/questions/11036689/solving-the-integer-knapsack
 * @param MaxWeight (in) maximum weight allowed in the knapsack
 * @param NumItems (in) number of items to select from for 
 *   inclusion in the knapsack
 * @param weight (in) pointer to array of size NumItems 
 *   with the weight of each item
 * @param value (in) pointer to array of size NumItems with the value 
 *   of each item
 * @return knapsack object
 */
Knapsack makeKnapsack(size_t MaxWeight, size_t NumItems, 
    size_t* weight, size_t* value);

/**
 * @brief destroy the knapsack object
 * @param k (in) knapsack object
 */
void destroyKnapsack(Knapsack k);

/**
 * @brief print the table used for computing the solution
 * @param k (in) knapsack object
 */
void printTable(Knapsack k);

/**
 * @brief get the items that form the solution
 * @remark runs in O(NumItems) time
 * @param k (in) knapsack object
 * @param size (inOut) length of solution array
 * @return array of item indexes (user must deallocate)
 */
size_t* getSolution(Knapsack k, size_t* size);

/**
 * @brief solve the knapsack problem
 * @remark runs in O(NumItems*MaxWeight) time
 * @param k (in) knapsack object
 * @return value of solution found
 */
size_t solve(Knapsack k);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif
