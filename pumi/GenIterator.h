/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef _GENERIC_ITERATOR_H_
#define _GENERIC_ITERATOR_H_

#ifdef __cplusplus

  template<class Iterator, class Entity> 
  class GenIterator
  {
  private:
    int dim;
    int topo; 
    Iterator theIter;
    Iterator iend; 
    Iterator ibegin;    
    Iterator prev; 
    void (*ptrFuncProcessingFilter) (Iterator&, Iterator&, void*, int, int); 
    void* ptr; 

  public:
    GenIterator(const Iterator &begin, 
	         const Iterator &end, 
		 int dim, 
		 int topo, 
		 void* ptr, 
		 void (*ptr2Func) (Iterator&, Iterator&, void* ptr, int, int)); 
    ~GenIterator() {}

    inline bool operator != (GenIterator const &) const;
    inline bool operator == (GenIterator const &) const;

   inline void* getPtr() const
   {
      return ptr; 
   }

   inline Entity * operator *() const
    {
      return *theIter; 
    }
    
    inline GenIterator & operator++() 
    {
      prev = theIter; 
      ++theIter;
      ptrFuncProcessingFilter(theIter, iend, ptr, dim, topo); 
      return *this;
    } 
    
    inline GenIterator operator++(int) 
    { 
      prev = theIter; 
      GenIterator tmp = *this; 
      ++(*this); 
      ptrFuncProcessingFilter(theIter, iend, ptr, dim, topo);
      return tmp; 
    }

     inline  bool end() const
     {
       return iend == theIter; 
     }

    inline void next ()
    {
       prev = theIter; 
       ++theIter; 
       ptrFuncProcessingFilter(theIter, iend, ptr, dim, topo); 
    }
   
    inline void reset ()
     {
         theIter = ibegin; 
	 prev = ibegin; 
     }

     inline void valid_check ()
     {
       if(theIter==iend)
       {
	 Iterator tmp = prev; 
	 ++tmp; 
	 ptrFuncProcessingFilter(tmp, iend, ptr, dim, topo);
	 if(tmp!=theIter)
	    theIter = tmp; 
       }
     }

    template<typename OtherIterator>
    inline const GenIterator& operator=(const GenIterator<OtherIterator, Entity>& other)
    {
      if(this==&other)
	 return *this; 

      dim = other.dim; 
      ptr = other.ptr; 
      theIter = other.theIter; 
      iend    = other.iend; 
      ibegin  = other.ibegin; 
      ptrFuncProcessingFilter = other.ptrFuncProcessingFilter; 
      prev = other.prev; 
      return *this; 
    }

    template<typename OtherIterator>
    GenIterator(const GenIterator<OtherIterator, Entity>& other)
    {
      dim = other.dim;
      ptr = other.ptr;
      theIter = other.theIter;
      iend = other.iend;
      ibegin = other.ibegin;
      prev = other.prev; 
      ptrFuncProcessingFilter = other.ptrFuncProcessingFilter; 
    }

};
 

template<typename Iterator, typename Entity>
inline bool GenIterator<Iterator, Entity> :: operator != (GenIterator const & other) const
  {
    return theIter != other.theIter;
  } 
  
  template<class Iterator, typename Entity>
  bool GenIterator<Iterator, Entity> :: operator == (GenIterator<Iterator, Entity> const & other) const
  {
    return theIter == other.theIter;
  } 

template<typename Iterator, typename Entity>
inline GenIterator<Iterator, Entity>::GenIterator(const Iterator &begin, const Iterator &end, int dim, int topo, void* argptr, void (*ptr2Func) (Iterator&, Iterator&, void*, int, int))
   :dim(dim), topo(topo), theIter(begin), iend(end), ibegin(begin), ptr(argptr)
{
   ptrFuncProcessingFilter = ptr2Func; 
   ptrFuncProcessingFilter(ibegin, iend, ptr, dim, topo);
   theIter = ibegin; 
   prev = ibegin; 
}

#endif
#endif
