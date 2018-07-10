#include "apfAggregateNumbering.h"
#include "apfAggregateNumberingClass.h"
#include "apfField.h"
#include "apfNumberingClass.h"
#include "apfTagData.h"
#include "apfUserData.h"
#include <PCU.h>
#include <pcu_util.h>
#include <set>
namespace apf
{
  // explicit instantiations
  template class AggregateNumberingOf<int>;
  template class AggregateNumberingOf<long>;
  // this assumes that the ranks returned by a sharing
  //  do not change when the pcu comm changes
  //  if the sharing can handle the pcu comm changing
  //  then this isn't needed
  class StaticSharing : public NormalSharing
  {
  protected:
    Sharing * shr;
    // original comm
    MPI_Comm old_cm;
    int old_rnk;
    // new comm
    MPI_Comm new_cm;
  public:
    StaticSharing(Mesh * m,
                  Sharing * s,
                  MPI_Comm ncm,
                  MPI_Comm pcm = MPI_COMM_NULL)
      : NormalSharing(m)
      , shr(s)
      , old_cm(pcm)
      , old_rnk(-1)
      , new_cm(ncm)
    {
      if(old_cm == MPI_COMM_NULL)
        old_cm = PCU_Get_Comm();
      new_cm = PCU_Get_Comm();
      old_rnk = PCU_Comm_Self();
    }
    virtual bool isOwned(MeshEntity * e)
    {
      int old_ownr = shr->getOwner(e);
      if(old_ownr == old_rnk)
        return true;
      return false;
    }
    virtual int getOwner(MeshEntity * e)
    {
      int old_ownr = shr->getOwner(e);
      return PCU_Foreign_To_Local(old_ownr,old_cm);
    }
    virtual void getCopies(MeshEntity * e, CopyArray & copies)
    {
      shr->getCopies(e,copies);
      for(int ii = 0; ii < copies.getSize(); ++ii)
        copies[ii].peer = PCU_Foreign_To_Local(copies[ii].peer,old_cm);
    }
    virtual bool isShared(MeshEntity * e)
    {
      return shr->isShared(e);
    }
  };
  template <class T>
  struct AggDataOfFunc : public FunctionBase<T>
  {
    AggDataOfFunc(AggregateNumberingOf<T> * n,
                  TagDataOf<T> * i)
      : num(n)
      , nd_ids(i)
    { }
    virtual void eval(MeshEntity * e, T * dat)
    {
      int nds = num->getShape()->countNodesOn(num->getMesh()->getType(e));
      int blks_per_nd = num->countBlocks();
      int dofs_per_blk = num->blockSize();
      int cmps = blks_per_nd * dofs_per_blk;
      int ownr = num->getSharing()->getOwner(e);
      int strd = num->getStride(ownr);
      int frst = num->getFirstDof(ownr);
      for(int nd = 0; nd < nds; ++nd)
      {
        T nd_id = -1;
        nd_ids->get(e,&nd_id);
        for(int blk = 0; blk < blks_per_nd; ++blk)
          for(int dof = 0; dof < dofs_per_blk; ++dof)
            dat[nd * cmps + blk * dofs_per_blk + dof] =
              frst + (nd_id * dofs_per_blk) + (blk * strd) + dof;
      }
    }
    AggregateNumberingOf<T> * num;
    TagDataOf<T> * nd_ids;
  };
  template <class T>
  T exscanT(T v);
  template <>
  int exscanT(int v)
  {
    return PCU_Exscan_Int(v);
  }
  template <>
  long exscanT(long v)
  {
    return PCU_Exscan_Long(v);
  }
  template <class T>
  void allgatherT(T i, T* o);
  template <>
  void allgatherT(int i, int* o)
  {
    PCU_Allgather_Int(i,o);
  }
  template <>
  void allgatherT(long i, long* o)
  {
    PCU_Allgather_Long(i,o);
  }
  template <class T>
  void AggregateNumberingOf<T>::init(const char * n,
                                     Mesh * m,
                                     Sharing * sh,
                                     FieldShape * s,
                                     int blks,
                                     int bs)
  {
    // trick AggFieldDataOf into only allocating one int per node
    shr = sh;
    //FieldBase::init(n,m,s,new AggDataOf<T>(this));
    NumberingOf<T>::components = 1;
    FieldBase::shape = s;
    FieldBase::mesh = m;
    FieldBase::name = n;
    nd_ids = new TagDataOf<T>;
    nd_ids->init(this);
    FieldBase::init(n,m,s,new UserDataBase<T>(new AggDataOfFunc<T>(this,static_cast<TagDataOf<T>*>(nd_ids))));
    NumberingOf<T>::components = blks * bs;
    blks_per_nd = blks;
    dofs_per_blk = bs;
    int dim = m->getDimension();
    for(int dd = 0; dd < dim; ++dd)
      lcl_strd += countOwnedNodes(m,dd,s,sh);
    int sz = PCU_Comm_Peers();
    slf = PCU_Comm_Self();
    // do we bother setting strds and frsts since they
    //  are reset in globalize which must be called after
    //  init for the numbering to be valid
    strds.resize(static_cast<unsigned>(sz));
    // all peers in the aggregation scope necessarily
    //  have the same stride
    int strd = PCU_Add_Int(lcl_strd);
    for(int ii = 0; ii < sz; ++ii)
      strds[ii] = strd;
    // the first dof id of the first local node is the # of
    //  nodes before the first local node times the # of
    //  dofs in the first block on a node
    T dof_id_offset = lcl_strd * dofs_per_blk;
    T frst = exscanT<T>(dof_id_offset);
    frsts.resize(static_cast<unsigned>(sz));
    allgatherT(frst,&frsts[0]);
    T nd_id = 0;
    for(int dd = 0; dd < dim; ++dd)
    {
      MeshIterator * it = m->begin(dd);
      MeshEntity * e = NULL;
      while((e = m->iterate(it)))
      {
        int nds = FieldBase::shape->countNodesOn(m->getType(e));
        // only really works if nds == 1 or 0
        //  due to the above comment ^
        for(int nd = 0; nd < nds; ++nd)
        {
          if(shr->isOwned(e))
          {
            nd_ids->set(e,&nd_id);
            ++nd_id;
          }
        }
      }
      m->end(it);
    }
  }
  template <class T>
  void AggregateNumberingOf<T>::init(Field * fld,
                                     Sharing * sh,
                                     int blks,
                                     int bs)
  {
    PCU_ALWAYS_ASSERT(fld->countComponents() == blks * bs);
    std::string nm = fld->getName();
    nm += "_num";
    init(nm.c_str(),
         fld->getMesh(),
         sh,
         fld->getShape(),
         blks,
         bs);
  }
  template <class T>
  void AggregateNumberingOf<T>::globalize()
  {
    synchronizeFieldData(nd_ids,shr,false);
    // any rank with frst == 0 was the zero rank
    //  in the aggregation scope under which init
    //  was called, the offsets from these
    //  zero-ranked processes are what we need to
    //  sum since they include in their strides
    //  all the nodes from the local aggregation scope
    int strd = getScopeStride();
    int frst = getLocalFirstDof();
    int scp_offset = (int)(!(bool)frst) * (strd * dofs_per_blk * blks_per_nd);
    int lcl_offset = PCU_Exscan_Int(scp_offset);
    lcl_offset += scp_offset; // make the excan a scan
    // remove the contribution from the local
    //  aggregation scope, because the scope can be
    //  any comm and only rank 0 in each scope
    //  contributes to the above scan
    lcl_offset -= (strd * dofs_per_blk * blks_per_nd);
    frst += lcl_offset;
    // strds and frsts needs to be resized/updated
    //  since we can have adjacent parts that were
    //  not in the aggregation scope under which
    //  init was called
    int sz = PCU_Comm_Peers();
    frsts.setSize(sz);
    allgatherT<T>(frst,&frsts[0]);
    // each aggregation scope has its own stride
    //  we store per-rank since that is vastly
    //  less complicated than keeping track
    //  of the initial scopes
    strds.setSize(sz);
    allgatherT<T>(strd,&strds[0]);
  }
  template <class T>
  void AggregateNumberingOf<T>::getAll(MeshEntity * e, T * dat)
  {
    reinterpret_cast<FieldDataOf<T>*>(FieldBase::data)->get(e,dat);
  }
  template <class T>
  T AggregateNumberingOf<T>::get(MeshEntity * e, int nd, int cmp)
  {
    int cmps = blks_per_nd * dofs_per_blk;
    int nds = FieldBase::shape->countNodesOn(FieldBase::mesh->getType(e));
    NewArray<T> data(nds*cmps);
    getAll(e,&data[0]);
    return data[nd*cmps + cmp];
  }
  AggNumbering * createAggNumbering(Field * f,
                                    int blocks,
                                    int dofs_per_block,
                                    MPI_Comm cm,
                                    Sharing * share)
  {
    bool dlt = false;
    if(!share)
    {
      share = getSharing(getMesh(f));
      dlt = true;
    }
    AggNumbering * n = new AggNumbering;
    MPI_Comm pcu_cm = PCU_Get_Comm();
    bool swtch = pcu_cm != cm;
    Sharing * init_share = share;
    if(swtch)
    {
      init_share = new StaticSharing(getMesh(f),share,cm);
      PCU_Switch_Comm(cm);
    }
    n->init(f,init_share,blocks,dofs_per_block);
    if(swtch)
    {
      PCU_Switch_Comm(pcu_cm);
      delete init_share;
    }
    // during init the static sharing may be used which
    //  only provides binary ownership information, not
    //  the owner rank of unowned entities, so we replace
    //  it with the passed-in sharing here before globalization
    //  which requires the rank information
    n->setSharing(share);
    n->globalize();
    f->getMesh()->addNumbering(n);
    return n;
  }
  AggNumbering * createAggNumbering(Mesh * m,
                                    const char * name,
                                    FieldShape * shape,
                                    int blocks,
                                    int dofs_per_block,
                                    MPI_Comm cm,
                                    Sharing * share)
  {
    bool dlt = false;
    if(!share)
    {
      share = getSharing(m);
      dlt = true;
    }
    AggNumbering * n = new AggNumbering;
    MPI_Comm pcu_cm = PCU_Get_Comm();
    bool swtch = pcu_cm != cm;
    Sharing * init_share = share;
    if(swtch)
    {
      init_share = new StaticSharing(m,share,cm);
      PCU_Switch_Comm(cm);
    }
    n->init(name,m,init_share,shape,blocks,dofs_per_block);
    if(swtch)
    {
      PCU_Switch_Comm(pcu_cm);
      delete init_share;
    }
    n->setSharing(share);
    n->globalize();
    m->addNumbering(n);
    return n;
  }
  /* Public API */
  Numbering * getNumbering(AggNumbering * n)
  {
    return static_cast<Numbering*>(n);
  }
  int countBlocks(AggNumbering * n)
  {
    return n->blockSize();
  }
  int countDOFsPerBlock(AggNumbering * n)
  {
    return n->blockSize();
  }
}
