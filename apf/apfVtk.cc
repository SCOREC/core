/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include <lionBase64.h>
#include <lionCompress.h>
#include "apfMesh.h"
#include "apfNumbering.h"
#include "apfNumberingClass.h"
#include "apfShape.h"
#include "apfFieldData.h"
#include <sstream>
#include <fstream>
#include <cassert>
#include <cstdlib>

namespace apf {

class HasAll : public FieldOp
{
  public:
    virtual bool inEntity(MeshEntity* e)
    {
      if (!f->getData()->hasEntity(e))
        ok = false;
      return false;
    }
    bool run(FieldBase* f_)
    {
      f = f_;
      ok = true;
      this->apply(f);
      return ok;
    }
  private:
    bool ok;
    FieldBase* f;
};

static bool isPrintable(FieldBase* f)
{
  HasAll op;
  return op.run(f);
}

static bool isNodal(FieldBase* f)
{
  return f->getShape() == f->getMesh()->getShape();
}

static bool isIP(FieldBase* f)
{
  FieldShape* s = f->getShape();
  Mesh* m = f->getMesh();
  int d = m->getDimension();
  for (int i=0; i < d; ++i)
    if (s->hasNodesIn(i))
      return false;
  return s->hasNodesIn(d);
}

static void describeArray(
    std::ostream& file,
    const char* name,
    int type,
    int size,
    bool isWritingBinary = false)
{
  file << "type=\"";
  const char* typeNames[3] = {"Float64","Int32","Int64"};
  file << typeNames[type];
  file << "\" Name=\"" << name;
  file << "\" NumberOfComponents=\"" << size;
  if (isWritingBinary)
  {
    file << "\" format=\"binary\"";
  }
  else
  {
    file << "\" format=\"ascii\"";  
  }
}

static void writePDataArray(
    std::ostream& file,
    const char* name,
    int type,
    int size,
    bool isWritingBinary = false)
{
  file << "<PDataArray ";
  describeArray(file,name,type,size,isWritingBinary);
  file << "/>\n";
}

static void writePDataArray(
    std::ostream& file,
    FieldBase* f,
    bool isWritingBinary = false)
{
  file << "<PDataArray ";
  describeArray(file,
      f->getName(),
      f->getScalarType(),
      f->countComponents(),
      isWritingBinary);
  file << "/>\n";
}

static void writePPoints(std::ostream& file,
    Field* f,
    bool isWritingBinary = false)
{
  file << "<PPoints>\n";
  writePDataArray(file,f,isWritingBinary);
  file << "</PPoints>\n";
}

static void writePPointData(std::ostream& file,
    Mesh* m, bool
    isWritingBinary = false)
{
  file << "<PPointData>\n";
  for (int i=0; i < m->countFields(); ++i)
  {
    Field* f = m->getField(i);
    if (isNodal(f) && isPrintable(f))
    {
      writePDataArray(file,f,isWritingBinary);
    }
  }
  for (int i=0; i < m->countNumberings(); ++i)
  {
    Numbering* n = m->getNumbering(i);
    if (isNodal(n) && isPrintable(n))
    {
      writePDataArray(file,n,isWritingBinary);
    }
  }
  for (int i=0; i < m->countGlobalNumberings(); ++i)
  {
    GlobalNumbering* n = m->getGlobalNumbering(i);
    if (isNodal(n) && isPrintable(n))
    {
      writePDataArray(file,n,isWritingBinary);
    }
  }
  file << "</PPointData>\n";
}

static int countIPs(FieldBase* f)
{
  /* assuming a non-mixed mesh here. for the already strained capabilities
     of VTK to accept IP fields, this is the best we can do */
  Mesh* m = f->getMesh();
  MeshIterator* it = m->begin(m->getDimension());
  MeshEntity* e = m->iterate(it);
  m->end(it);
  return f->countNodesOn(e);
}

static std::string getIPName(FieldBase* f, int point)
{
  std::stringstream ss;
  /* People looking at these files get scared of 0-based indexing
                                     V                        */
  ss << f->getName() << '_' << (point+1);
  return ss.str();
}

static void writeIP_PCellData(std::ostream& file,
    FieldBase* f,
    bool isWritingBinary = false)
{
  int n = countIPs(f);
  for (int p=0; p < n; ++p)
  {
    std::string s= getIPName(f,p);
    writePDataArray(file,
        s.c_str(),f->getScalarType(),
        f->countComponents(),
        isWritingBinary);
  }
}

static void writePCellParts(std::ostream& file, bool isWritingBinary = false)
{
  writePDataArray(file, "apf_part", apf::Mesh::INT, 1, isWritingBinary);
}

static void writePCellData(std::ostream& file,
    Mesh* m,
    bool isWritingBinary = false)
{
  file << "<PCellData>\n";
  for (int i=0; i < m->countFields(); ++i)
  {
    Field* f = m->getField(i);
    if (isIP(f) && isPrintable(f))
    {
      writeIP_PCellData(file,f,isWritingBinary);
    }
  }
  for (int i=0; i < m->countNumberings(); ++i)
  {
    Numbering* n = m->getNumbering(i);
    if (isIP(n) && isPrintable(n))
    {
      writeIP_PCellData(file,n,isWritingBinary);
    }
  }
  for (int i=0; i < m->countGlobalNumberings(); ++i)
  {
    GlobalNumbering* n = m->getGlobalNumbering(i);
    if (isIP(n) && isPrintable(n))
    {
      writeIP_PCellData(file,n,isWritingBinary);
    }
  }
  writePCellParts(file,isWritingBinary);
  file << "</PCellData>\n";
}

static std::string getPieceFileName(const char* prefix, int id)
{
  std::stringstream ss;
  ss << prefix << id << ".vtu";
  return ss.str();
}

static std::string stripPath(std::string const& s)
{
  size_t i = s.rfind('/');
  if (i == std::string::npos)
  {
    return s;
  }
  return s.substr(i + 1, std::string::npos);
}

static void writePSources(std::ostream& file, const char* prefix)
{
  for (int i=0; i < PCU_Comm_Peers(); ++i)
  {
    std::string fileName = stripPath(getPieceFileName(prefix,i));
    file << "<Piece Source=\"" << fileName << "\"/>\n";
  }
}

static void writePvtuFile(const char* prefix,
    Mesh* m,
    bool isWritingBinary = false)
{
  std::string fileName = prefix;
  fileName += ".pvtu";
  std::ofstream file(fileName.c_str());
  assert(file.is_open());
  file << "<VTKFile type=\"PUnstructuredGrid\">\n";
  file << "<PUnstructuredGrid GhostLevel=\"0\">\n";
  writePPoints(file,m->getCoordinateField(),isWritingBinary);
  writePPointData(file,m,isWritingBinary);
  writePCellData(file,m,isWritingBinary);
  writePSources(file,prefix);
  file << "</PUnstructuredGrid>\n";
  file << "</VTKFile>\n";
}

static void writeDataHeader(std::ostream& file,
    const char* name,
    int type,
    int size,
    bool isWritingBinary = false)
{
  file << "<DataArray ";
  describeArray(file,name,type,size,isWritingBinary);
  file << ">\n";
}

static void writeEncodedArray(std::ostream& file,
    unsigned int dataLenBytes,
    char* dataToEncode)
{
  if ( lion::can_compress )
  {
    //build the data header and compress dataToEncode
    long lensToEncode[4];
    lensToEncode[0] = 1; //data is compressed in one block
    lensToEncode[1] = dataLenBytes; //size of each block before compression
    lensToEncode[2] = dataLenBytes; //size of the final block before compression
    unsigned long dataCompressedLen = lion::compressBound(dataLenBytes);
        //gets the size to allocate for the compressed data array
    char* dataCompressed = new char[dataCompressedLen]();
    lion::compress( (void*)dataCompressed,
        dataCompressedLen,
        (void*)dataToEncode,
        dataLenBytes);
        //compresses and sets dataCompressedLen to the size of the correct length
    lensToEncode[3] = dataCompressedLen; //size of the compressed block
    //encode and output the compressed data
    file << lion::base64Encode( (char*)lensToEncode, sizeof(lensToEncode) );
    file << lion::base64Encode( dataCompressed, dataCompressedLen ) << '\n';
    delete[] dataCompressed;
  }
  else
  {
    //not compressing, encode and output
    file << lion::base64Encode( (char*)&dataLenBytes, sizeof(dataLenBytes) );
    file << lion::base64Encode( dataToEncode, dataLenBytes ) << '\n';
  }
}

template <class T>
static void writeNodalField(std::ostream& file,
    FieldBase* f,
    DynamicArray<Node>& nodes,
    bool isWritingBinary = false)
{
  int nc = f->countComponents();
  writeDataHeader(file,f->getName(),f->getScalarType(),nc,isWritingBinary);
  NewArray<T> nodalData(nc);
  FieldDataOf<T>* data = static_cast<FieldDataOf<T>*>(f->getData());
  if (isWritingBinary)
  {
    unsigned int dataLen = nc * nodes.getSize();
    unsigned int dataLenBytes = dataLen*sizeof(T);
    T* dataToEncode = new T[dataLen]();
    unsigned int dataIndex = 0;
    for (size_t i = 0; i < nodes.getSize(); ++i)
    {
      data->getNodeComponents(nodes[i].entity,nodes[i].node,&(nodalData[0]));
      for (int j = 0; j < nc; ++j)
      {
        dataToEncode[dataIndex] = nodalData[j];
        dataIndex++;
      }
    }
    writeEncodedArray(file, dataLenBytes, (char*)dataToEncode);
    delete [] dataToEncode;
  }
  else
  {
    for (size_t i = 0; i < nodes.getSize(); ++i)
    {
      data->getNodeComponents(nodes[i].entity,nodes[i].node,&(nodalData[0]));
      for (int j = 0; j < nc; ++j)
      {
        file << nodalData[j] << ' ';
      }
      file << '\n';
    }
  }
  file << "</DataArray>\n"; 
}

static void writePoints(std::ostream& file,
    Mesh* m,
    DynamicArray<Node>& nodes,
    bool isWritingBinary = false)
{
  file << "<Points>\n";
  writeNodalField<double>(file,m->getCoordinateField(),nodes,isWritingBinary);
  file << "</Points>\n"; 
}

static int countElementNodes(Numbering* n, MeshEntity* e)
{
  return n->getShape()->getEntityShape(n->getMesh()->getType(e))->countNodes();
}

static void writeConnectivity(std::ostream& file,
    Numbering* n,
    bool isWritingBinary = false)
{
  file << "<DataArray type=\"Int32\" Name=\"connectivity\"";
  if (isWritingBinary)
  {
    file << " format=\"binary\"";
  }
  else
  {
    file << " format=\"ascii\"";
  }
  file << ">\n";
  Mesh* m = n->getMesh();
  MeshEntity* e;
  if (isWritingBinary)
  {
    // TODO: see if we can do this with only one loop
    MeshIterator* elements = m->begin(m->getDimension());
    unsigned int dataLen = 0;
    while ((e = m->iterate(elements)))
    {
      dataLen += countElementNodes(n,e);
    }
    m->end(elements);
    unsigned int dataLenBytes = dataLen*sizeof(int);
    int* dataToEncode = new int[dataLen]();
    elements = m->begin(m->getDimension());
    unsigned int dataIndex = 0;
    while ((e = m->iterate(elements)))
    {
      int nen = countElementNodes(n,e);
      NewArray<int> numbers(nen);
      getElementNumbers(n,e,numbers);
      for (int i=0; i < nen; ++i)
      {
        dataToEncode[dataIndex] = numbers[i];
        dataIndex++;
      }
    }
    m->end(elements);
    writeEncodedArray(file, dataLenBytes, (char*)dataToEncode);
    delete [] dataToEncode;
  }
  else
  {
    MeshIterator* elements = m->begin(m->getDimension());
    while ((e = m->iterate(elements)))
    {
      int nen = countElementNodes(n,e);
      NewArray<int> numbers(nen);
      getElementNumbers(n,e,numbers);
      for (int i=0; i < nen; ++i)
      {
        file << numbers[i] << ' ';
      }
      file << '\n';
    }
    m->end(elements);
  }
  file << "</DataArray>\n";
}

static void writeOffsets(std::ostream& file,
    Numbering* n,
    bool isWritingBinary = false)
{
  file << "<DataArray type=\"Int32\" Name=\"offsets\"";
  if (isWritingBinary)
  {
    file << " format=\"binary\"";
  }
  else
  {
    file << " format=\"ascii\"";
  }
  file << ">\n";
  Mesh* m = n->getMesh();
  MeshEntity* e;
  if (isWritingBinary)
  {
    // TODO: see if we can do this with only one loop
    MeshIterator* elements = m->begin(m->getDimension());
    unsigned int dataLen = 0;
    while ((e = m->iterate(elements)))
    {
      dataLen++;
    }
    m->end(elements);
    unsigned int dataLenBytes = dataLen*sizeof(int);
    int* dataToEncode = new int[dataLen]();
    elements = m->begin(m->getDimension());
    unsigned int dataIndex = 0;
    int offset = 0;
    while ((e = m->iterate(elements)))
    {
      offset += countElementNodes(n,e);
      dataToEncode[dataIndex] = offset;
      dataIndex++;
    }
    m->end(elements);
    writeEncodedArray(file, dataLenBytes, (char*)dataToEncode);
    delete [] dataToEncode;
  }
  else
  {
    MeshIterator* elements = m->begin(m->getDimension());
    int offset = 0;
    while ((e = m->iterate(elements)))
    {
      offset += countElementNodes(n,e);
      file << offset << '\n';
    }
    m->end(elements);
  }
  file << "</DataArray>\n";
}

static void writeTypes(std::ostream& file,
    Mesh* m,
    bool isWritingBinary = false)
{
  file << "<DataArray type=\"Int32\" Name=\"types\"";
  if (isWritingBinary)
  {
    file << " format=\"binary\"";
  }
  else
  {
    file << " format=\"ascii\"";
  }
  file << ">\n";
  MeshEntity* e;
  int order = m->getShape()->getOrder();
  static int vtkTypes[Mesh::TYPES][2] =
    /* order
       linear,quadratic
       V  V */
    {{ 1,-1}//vertex
    ,{ 3,21}//edge
    ,{ 5,22}//triangle
    ,{ 9,23}//quad
    ,{10,24}//tet
    ,{12,25}//hex
    ,{13,-1}//prism
    ,{14,-1}//pyramid
  };
  if (isWritingBinary)
  {
    unsigned int dataLen = 0;
    MeshIterator* elements = m->begin(m->getDimension());
    while ((e = m->iterate(elements)))
    {
      dataLen++;
    }
    m->end(elements);
    unsigned int dataLenBytes = dataLen*sizeof(int);
    int* dataToEncode = new int[dataLen]();
    elements = m->begin(m->getDimension());
    unsigned int dataIndex = 0;
    while ((e = m->iterate(elements)))
    {
      dataToEncode[dataIndex] = vtkTypes[m->getType(e)][order-1];
      dataIndex++;
    }
    m->end(elements);
    writeEncodedArray(file, dataLenBytes, (char*)dataToEncode);
    delete [] dataToEncode;
  }
  else
  {
    MeshIterator* elements = m->begin(m->getDimension());
    while ((e = m->iterate(elements)))
    {
      file << vtkTypes[m->getType(e)][order-1] << '\n';
    }
    m->end(elements);
  }
  file << "</DataArray>\n";
}

static void writeCells(std::ostream& file,
    Numbering* n,
    bool isWritingBinary = false)
{
  file << "<Cells>\n";
  writeConnectivity(file,n,isWritingBinary);
  writeOffsets(file,n,isWritingBinary);
  writeTypes(file,n->getMesh(),isWritingBinary);
  file << "</Cells>\n"; 
}

static void writePointData(std::ostream& file,
    Mesh* m,
    DynamicArray<Node>& nodes,
    bool isWritingBinary = false)
{
  file << "<PointData>\n";
  for (int i=0; i < m->countFields(); ++i)
  {
    Field* f = m->getField(i);
    if (isNodal(f) && isPrintable(f))
    {
      writeNodalField<double>(file,f,nodes,isWritingBinary);
    }
  }
  for (int i=0; i < m->countNumberings(); ++i)
  {
    Numbering* n = m->getNumbering(i);
    if (isNodal(n) && isPrintable(n))
    {
      writeNodalField<int>(file,n,nodes,isWritingBinary);
    }
  }
  for (int i=0; i < m->countGlobalNumberings(); ++i)
  {
    GlobalNumbering* n = m->getGlobalNumbering(i);
    if (isNodal(n) && isPrintable(n))
    {
      writeNodalField<long>(file,n,nodes,isWritingBinary);
    }
  }
  file << "</PointData>\n";
}

template <class T>
class WriteIPField : public FieldOp
//TODO change to write encoded data
{
  public:
    int point;
    int components;
    NewArray<T> ipData;
    FieldDataOf<T>* data;
    MeshEntity* entity;
    std::ostream* fp;
    virtual bool inEntity(MeshEntity* e)
    {
      entity = e;
      return true;
    }
    virtual void atNode(int node)
    {
      if (node != point)
        return;
      data->getNodeComponents(entity,node,&(ipData[0]));
      for (int i=0; i < components; ++i)
        (*fp) << ipData[i] << ' ';
      (*fp) << '\n';
    }
    void runOnce(FieldBase* f)
    {
      std::string s = getIPName(f,point);
      components = f->countComponents();
      writeDataHeader(*fp,s.c_str(),f->getScalarType(),f->countComponents());
      ipData.allocate(components);
      data = static_cast<FieldDataOf<T>*>(f->getData());
      apply(f);
      (*fp) << "</DataArray>\n"; 
    }
    void run(std::ostream& file, FieldBase* f)
    {
      fp = &file;
      int n = countIPs(f);
      for (point=0; point < n; ++point)
        runOnce(f);
    }
};

static void writeCellParts(std::ostream& file, 
    Mesh* m,
    bool isWritingBinary = false)
{
  writeDataHeader(file, "apf_part", apf::Mesh::INT, 1, isWritingBinary);
  size_t n = m->count(m->getDimension());
  int id = m->getId();
  if (isWritingBinary)
  {
    unsigned int dataLen = n;
    unsigned int dataLenBytes = dataLen*sizeof(int);
    int* dataToEncode = new int[dataLen]();
    for (size_t i = 0; i < n; ++i )
    {
      dataToEncode[i] = id;
    }
    writeEncodedArray(file, dataLenBytes, (char*)dataToEncode);
    file << "</DataArray>\n";
    delete [] dataToEncode;
  }
  else
  { 
    for (size_t i = 0; i < n; ++i)
    {
      file << id << '\n';
    }
    file << "</DataArray>\n";
  }
}

static void writeCellData(std::ostream& file,
    Mesh* m,
    bool isWritingBinary = false)
{
  file << "<CellData>\n";
  WriteIPField<double> wd;
  for (int i=0; i < m->countFields(); ++i)
  {
    Field* f = m->getField(i);
    if (isIP(f) && isPrintable(f))
    {
      wd.run(file,f);
    }
  }
  WriteIPField<int> wi;
  for (int i=0; i < m->countNumberings(); ++i)
  {
    Numbering* n = m->getNumbering(i);
    if (isIP(n) && isPrintable(n))
    {
      wi.run(file,n);
    }
  }
  WriteIPField<long> wl;
  for (int i=0; i < m->countGlobalNumberings(); ++i)
  {
    GlobalNumbering* n = m->getGlobalNumbering(i);
    if (isIP(n) && isPrintable(n))
    {
      wl.run(file,n);
    }
  }
  writeCellParts(file, m, isWritingBinary);
  file << "</CellData>\n";
}

//function checks if the machine is a big or little endian machine
bool isBigEndian()
{
  union {
      unsigned int i;
      char c[4];
  } bint = {0x01020304};
  return bint.c[0] == 1;
}

static void writeVtuFile(const char* prefix,
    Numbering* n,
    bool isWritingBinary = false)
{
  PCU_Barrier();
  double t0 = PCU_Time();
  std::string fileName = getPieceFileName(prefix,PCU_Comm_Self());
  std::stringstream buf;
  Mesh* m = n->getMesh();
  DynamicArray<Node> nodes;
  getNodes(n,nodes);
  buf << "<VTKFile type=\"UnstructuredGrid\"";
  if (isWritingBinary)
  {
    buf << " byte_order=";
    if (isBigEndian())
    {
      buf << "\"BigEndian\"";
    }
    else
    {
      buf << "\"LittleEndian\"";
    }
    if (lion::can_compress )
    {
      //TODO determine what the header_type should be definatively
      buf << " header_type=\"UInt64\"";
      buf << " compressor=\"vtkZLibDataCompressor\"";
    }
    else
    {
      buf << " header_type=\"UInt32\"";
    }
  }
  buf<< ">\n";
  buf << "<UnstructuredGrid>\n";
  buf << "<Piece NumberOfPoints=\"" << nodes.getSize();
  buf << "\" NumberOfCells=\"" << m->count(m->getDimension());
  buf << "\">\n";
  writePoints(buf,m,nodes,isWritingBinary);
  writeCells(buf,n,isWritingBinary);
  writePointData(buf,m,nodes,isWritingBinary);
  writeCellData(buf,m,isWritingBinary);
  buf << "</Piece>\n";
  buf << "</UnstructuredGrid>\n";
  buf << "</VTKFile>\n";
  PCU_Barrier();
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    printf("writeVtuFile into buffers: %f seconds\n", t1 - t0);
  }
  { //block forces std::ofstream destructor call
    std::ofstream file(fileName.c_str());
    assert(file.is_open());
    file << buf.rdbuf();
  }
  PCU_Barrier();
  double t2 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    printf("writeVtuFile buffers to disk: %f seconds\n", t2 - t1);
  }
}

void writeVtkFiles(const char* prefix, Mesh* m)
{
//*** this function is now writing ASCII encoded, uncompressed vtk files
  writeASCIIVtkFiles(prefix, m);
  // --------------------------ignore for now----------------------------------
  //*** this function writes vtk files with binary encoding ***
  //use writeASCIIVtkFiles for ASCII encoding (not recommended)
  // bool isWritingBinary = true;
  // double t0 = PCU_Time();
  // if (!PCU_Comm_Self())
  // {
  //   writePvtuFile(prefix, m, isWritingBinary);
  // }
  // Numbering* n = numberOverlapNodes(m,"apf_vtk_number");
  // m->removeNumbering(n);
  // writeVtuFile(prefix, n, isWritingBinary);
  // double t1 = PCU_Time();
  // if (!PCU_Comm_Self())
  // {
  //   printf("vtk files %s written in %f seconds\n", prefix, t1 - t0);
  // }
  // delete n;
  // -------------------------------------------------------------------------
}

void writeOneVtkFile(const char* prefix, Mesh* m)
{
  /* creating a non-collective numbering is
     a tad bit risky, but we should be fine
     given the current state of the code */
  
  // bool isWritingBinary = true;
  Numbering* n = numberOverlapNodes(m,"apf_vtk_number");
  m->removeNumbering(n);
  // writeVtuFile(prefix, n, isWritingBinary);
  writeVtuFile(prefix, n);
  delete n;
}

void writeASCIIVtkFiles(const char* prefix, Mesh* m)
{
  //*** this function writes vtk files with ASCII encoding ***
  //*** not recommended, use writeVtkFiles instead ***
  double t0 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    writePvtuFile(prefix, m);
  }
  Numbering* n = numberOverlapNodes(m,"apf_vtk_number");
  m->removeNumbering(n);
  writeVtuFile(prefix, n);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    printf("vtk files %s written in %f seconds\n", prefix, t1 - t0);
  }
  delete n;
}

void writeBinaryVtkFiles(const char* prefix, Mesh* m)
{
  //*** this function writes vtk files with binary encoding ***
  //use writeASCIIVtkFiles for ASCII encoding (not recommended)
  bool isWritingBinary = true;
  double t0 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    writePvtuFile(prefix, m, isWritingBinary);
  }
  Numbering* n = numberOverlapNodes(m,"apf_vtk_number");
  m->removeNumbering(n);
  writeVtuFile(prefix, n, isWritingBinary);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
  {
    printf("vtk files %s written in %f seconds\n", prefix, t1 - t0);
  }
  delete n;
}

}
