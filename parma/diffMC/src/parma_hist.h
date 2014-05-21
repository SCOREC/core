#ifndef PARMA_HIST_H_
#define PARMA_HIST_H_

#include <vector>
#include <string>

class hist {
  public:
    hist();
    ~hist();
    void add(const int p);
    int getSz();
    void print(const int n, std::string outfile);
    void clear();
  private:
    std::vector<int> d;
};

#endif


