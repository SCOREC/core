#ifndef PARMA_PRIORITY_H_
#define PARMA_PRIORITY_H_

class priorityList {
   public:
      int entDim[4];
      int priority[4];  // not active if < 0
      void print();
      void sort(int (*inPrio)[4]);
};         

#endif
