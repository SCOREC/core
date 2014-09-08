#include "apfArrayData.h"
#include "apfNumbering.h"
#include "apfTagData.h"

namespace apf {

template <class T>
class ArrayDataOf : public FieldDataOf<T>
{
  public:
    virtual void init(FieldBase* f)
    {
      /* this class inherits a variable (field),
         lets initialize it */
      this->field = f;
      /* this has to set up the array */
      FieldShape* s = f->getShape();
      const char* name = s->getName();
      Numbering* n = f->getMesh()->findNumbering(name);
      if (n==NULL)
        num_var = numberOverlapNodes(f->getMesh(),name,s);
      else
        num_var = n;
      arraySize = f->countComponents()*countNodes(num_var);
      dataArray = new T[arraySize];
    }
    virtual ~ArrayDataOf()
    {
      /* this has to destroy the array */
      delete [] dataArray;
    }
    virtual bool hasEntity(MeshEntity*)
    {
      /* this should get a bit more complex:
         return true if the field has nodes on (e) */
      return true;
    }
    virtual void removeEntity(MeshEntity*)
    {
      /* this will remain an empty function...
         I don't think  we want to remove entities from frozen fields */
      fail("removeEntity called on frozen field data");
    }
    virtual void get(MeshEntity* e, T* data)
    {
      /* this retrieves all the data associated with (e) */
      int first_node_index = getNumber(this->num_var,e,0,0);
      int num_nodes = this->field->countNodesOn(e);
      int num_components = this->field->countComponents();
      int start = first_node_index*num_components;
      for (int i=0; i<num_nodes*num_components; i++) {
         data[i] = this->dataArray[start+i];
      }
    }
    virtual void set(MeshEntity* e, T const* data)
    {
      /* this stores all the data associated with (e) */
      int first_node_index = getNumber(this->num_var,e,0,0);
      int num_nodes = this->field->countNodesOn(e);
      int num_components = this->field->countComponents();
      int start = first_node_index*num_components;
      for (int i=0; i<num_nodes*num_components; i++) {
          this->dataArray[start+i] = data[i];
      }
    }

    virtual bool isFrozen() {
      return true;
    }

    T* getDataArray() {
      return this->dataArray;
    }

  private:
    /* data variables go here */
    Numbering* num_var; 
    int arraySize;
    T* dataArray;
};

template <class T>
void freezeFieldData(FieldBase* field)
{
  /* make a new data store of array type */
  ArrayDataOf<T>* newData = new ArrayDataOf<T>();
  /* call the init function to setup storage */
  newData->init(field);
  /* get the old data store */
  FieldDataOf<T>* oldData = static_cast<FieldDataOf<T>*>(field->getData());
  /* call the set function to fill with values */
  copyFieldData<T>(newData,oldData);
  /* replace the old data store with this one */
  field->changeData(newData);
}

template <class T>
void unfreezeFieldData(FieldBase* field) {
  // make a new data store of tag type
  TagDataOf<T>* newData = new TagDataOf<T>();
  // call init function to setup storage
  newData->init(field);
  // get the old data store
  FieldDataOf<T>* oldData = static_cast<FieldDataOf<T>*>(field->getData());
  // call set function to fill with values
  copyFieldData<T>(newData,oldData);
  // replace old data store with this one
  field->changeData(newData);
}

/* instantiate here */
template void freezeFieldData<int>(FieldBase* field);
template void freezeFieldData<double>(FieldBase* field);
template void unfreezeFieldData<int>(FieldBase* field);
template void unfreezeFieldData<double>(FieldBase* field);

double* getArrayData(Field* f) {
  if (!isFrozen(f)) {
    return 0;
  } else {
    FieldDataOf<double>* p = f->getData();
    ArrayDataOf<double>* a = static_cast<ArrayDataOf<double>* > (p);
    return a->getDataArray();
  }
}

}
