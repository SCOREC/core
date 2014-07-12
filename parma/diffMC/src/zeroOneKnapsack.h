#ifndef _KNAPSACK_H_
#define _KNAPSACK_H_

// http://stackoverflow.com/questions/11036689/solving-the-integer-knapsack

#include <iostream>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <math.h>

/**
 * @brief implements the zero-one knapsack problem
 * @remark code based on http://stackoverflow.com/questions/11036689/solving-the-integer-knapsack
 */
class knapsack {
public:
    /**
     * @brief create the knapsack object
     * @param MaxWeight (in) maximum weight allowed in the knapsack
     * @param NumItems (in) number of items to select from for inclusion in the knapsack
     * @param weight (in) pointer to integer array of size NumItems with the weight of each item
     * @param value (in) pointer to integer array of size NumItems with the value of each item
     */
    knapsack(int MaxWeight, int NumItems, int* weight, int* value);

    ~knapsack();

    /**
     * @brief print the table used for computing the solution
     * @return 0 on success, non-zero on failure
     */
    int printTable();

    /**
     * @brief get the items that form the solution
     * @remark runs in O(NumItems) time
     * @param solnIdx (InOut) vector of item indexes
     * @return 0 on success, non-zero on failure
     */
    int getSolution(std::vector<int>& solnIdx);

    /**
     * @brief solve the knapsack problem
     * @remark runs in O(NumItems*MaxWeight) time
     * @return value of solution found
     */
    int solve();
private:
    int **M;
    int *weight;
    int *value;
    int maxWeight;
    int numItems;
};

#endif
