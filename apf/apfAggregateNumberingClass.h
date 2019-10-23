#ifndef APF_AGGREGATE_NUMBERING_CLASS_H_
#define APF_AGGREGATE_NUMBERING_CLASS_H_
#include "apfFieldData.h"
#include "apfNumberingClass.h"
namespace apf
{
  template <class T>
  class AggregateNumberingOf : public NumberingOf<T>
  {
  public:
    AggregateNumberingOf()
      : NumberingOf<T>()
      , nd_ids(NULL)
      , shr(NULL)
      , blks_per_nd(0)
      , dofs_per_blk(0)
      , slf(-1)
      , strds()
      , frsts()
      , lcl_strd(0)
    { }
    void init(const char * n,
              Mesh * m,
              Sharing * sh,
              FieldShape * s,
              int blks,
              int bs);
    void init(Field * f, Sharing * sh, int blks, int bs);
    void globalize();
    virtual void getAll(MeshEntity * e, T * dat);
    virtual T get(MeshEntity * e, int nd, int cmp);
    virtual void set(MeshEntity*,int,int,T) {}
    Sharing * getSharing() const { return shr; }
    void setSharing(Sharing * s) { shr = s; }
    int countBlocks() const { return blks_per_nd; }
    int blockSize() const { return dofs_per_blk; }
    T getLocalStride() const { return lcl_strd; }
    T getScopeStride() const { return strds[slf]; }
    T getLocalFirstDof() const { return frsts[slf]; }
    T getStride(int peer) const { return strds[peer]; };
    T getFirstDof(int peer) const { return frsts[peer]; }
  private:
    FieldDataOf<T> * nd_ids;
    Sharing * shr;
    int blks_per_nd;
    int dofs_per_blk;
    int slf;
    DynamicArray<T> strds;
    DynamicArray<T> frsts;
    int lcl_strd;
  };
}
#endif
