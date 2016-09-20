/****************************************************************************** 

  (c) 2012-2016 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.
 
*******************************************************************************/
#ifndef PUMI_LIST_H
#define PUMI_LIST_H

#include <iterator>
#include <cstddef>

class ListMember
{
  public:
    ListMember():previous(this),next(this) {}
    void linkAfter(ListMember* p)
    {
      this->next = p->next;
      p->next = this;
      this->previous = p;
      this->next->previous = this;
    }
    void unlink()
    {
      this->previous->next = this->next;
      this->next->previous = this->previous;
    }
    ListMember* previous;
    ListMember* next;
};

template <class T>
class ListIterator : public std::iterator<std::forward_iterator_tag,T*>
{
  public:
    ListIterator():current(0) {}
    ListIterator(ListMember* p):current(p) {}
    bool operator==(ListIterator<T> const& other) const
    {
      return current == other.current;
    }
    bool operator!=(ListIterator<T> const& other) const
    {
      return current != other.current;
    }
    ListIterator<T>& operator++()
    {
      current = current->next;
      return *this;
    }
    ListIterator<T> operator++(int)
    {
      ListIterator<T> temp(*this);
      ++(*this);
      return temp;
    }
    T* operator*() const
    {
      return static_cast<T*>(current);
    }
  private:
    ListMember* current;
};

class List
{
  public:
    List():count(0) {}
    void push_back(ListMember* m)
    {
      m->linkAfter(head.previous);
      ++count;
    }
    void remove(ListMember* m)
    {
      m->unlink();
      --count;
    }
    template <class T>
    ListIterator<T> begin() {return ListIterator<T>(head.next);}
    template <class T>
    ListIterator<T> end() {return ListIterator<T>(&head);}
    std::size_t size() const {return count;}
  private:
    std::size_t count;
    ListMember head;
};


#endif
