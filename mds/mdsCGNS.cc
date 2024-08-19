/* Author: Dr Andrew Parker (2019) - FGE Ltd 

   Notes: these functions have been exclusively developed using meshes from Salome

*/

#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "apfShape.h"
#include "gmi.h" /* this is for gmi_getline... */
#include <lionPrint.h>
#include <PCU_C.h>
#include <apf.h>
#include <apfConvert.h>
#include "apfFieldData.h"
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfNumbering.h>

#include <cstdio>
#include <cstring>
#include <pcu_util.h>
#include <cstdlib>

#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
#include <list>
#include <algorithm>
#include <unordered_set>
#include <utility>
//
#ifdef HAVE_CGNS
//
#include <cgns_io.h>
#include <pcgnslib.h>
//
#endif

namespace
{
#ifndef NDEBUG                 // debug settings, cmake double negative....
const bool debugOutput = true; // probably will not get away with c++17
//static constexpr bool debugOutput = true; // probably will not get away with c++17
#else // optimised setting
const bool debugOutput = false; // probably will not get away with c++17
                                //static constexpr bool debugOutput = false; // probably will not get away with c++17
#endif

static std::string SupportedCGNSElementTypeToString(const CGNS_ENUMT(ElementType_t) & elementType)
{
  if (elementType == CGNS_ENUMV(NODE))
    return "NODE";
  else if (elementType == CGNS_ENUMV(BAR_2))
    return "BAR_2";
  else if (elementType == CGNS_ENUMV(QUAD_4))
    return "QUAD_4";
  else if (elementType == CGNS_ENUMV(TRI_3))
    return "TRI_3";
  else if (elementType == CGNS_ENUMV(HEXA_8))
    return "HEXA_8";
  else if (elementType == CGNS_ENUMV(TETRA_4))
    return "TETRA_4";
  else if (elementType == CGNS_ENUMV(PYRA_5))
    return "PYRA_5";
  else
  {
    std::cout << __LINE__ << " "
              << "Element type "
              << " " << cg_ElementTypeName(elementType) << " not supported by apf/PUMI " << std::endl;
    return "";
  }
}

struct MeshData
{
  MeshData(const int &s, const int &f, const CGNS_ENUMT(DataType_t) & dt, const CGNS_ENUMT(GridLocation_t) & l, const std::string &n) : si(s), fi(f), datatype(dt), location(l), name(n)
  {
  }

  int si = -1;
  int fi = -1;
  CGNS_ENUMT(DataType_t)
  datatype;
  CGNS_ENUMT(GridLocation_t)
  location;
  std::string name;
  bool process = true;

  bool operator==(const MeshData &other) const
  {
    if (si != other.si)
      return false;

    if (fi != other.fi)
      return false;

    if (datatype != other.datatype)
      return false;

    if (location != other.location)
      return false;

    if (name != other.name)
      return false;

    return true;
  }

  bool operator!=(const MeshData &other) const
  {
    return !this->operator==(other);
  }
};

struct MeshDataGroup
{
  std::map<int, MeshData> components;
  bool process = true;

  CGNS_ENUMT(GridLocation_t)
  location() const
  {
    return components.at(0).location;
  }

  CGNS_ENUMT(DataType_t)
  datatype() const
  {
    return components.at(0).datatype;
  }

  std::string name() const
  {
    // pattern matching component writer in apfCGNS.cc, patterns must be kept in-sync with this file
    const std::string end("_[");
    const std::size_t found = components.at(0).name.find(end);
    if (found != std::string::npos)
      return components.at(0).name.substr(0, found);
    else if (components.size() == 1)
      return components.at(0).name;
    else
    {
      info();
      PCU_ALWAYS_ASSERT_VERBOSE(true == false, "Problem with name");
      return "";
    }
  }

  std::size_t size() const
  {
    return components.size();
  }

  int sIndex(const int index) const
  {
    return components.at(index).si;
  }

  int fIndex(const int index) const
  {
    return components.at(index).fi;
  }

  void insert(int i, MeshData &d)
  {
    components.insert(std::make_pair(i, d));
  }

  void info() const
  {
    if (components.size() == 1)
    {
      std::cout << "Scalar Group has " << components.size() << " related componenets: " << std::endl;
      for (const auto m : components)
        std::cout << "Field " << m.second.name << " @ " << m.second.si << " " << m.second.fi << std::endl;
    }
    else if (components.size() == 3)
    {
      std::cout << "Vector Group has " << components.size() << " related componenets: " << std::endl;
      for (const auto m : components)
        std::cout << "Field " << m.second.name << " @ " << m.second.si << " " << m.second.fi << std::endl;
    }
    else if (components.size() == 9)
    {
      std::cout << "Matrix Group has " << components.size() << " related componenets: " << std::endl;
      for (const auto m : components)
        std::cout << "Field " << m.second.name << " @ " << m.second.si << " " << m.second.fi << std::endl;
    }
    else
    {
      PCU_ALWAYS_ASSERT_VERBOSE(true == false, "Tensor not accounted for");
    }
  }
};

template <typename Arg, typename... Args>
void DebugParallelPrinter(PCU_t h, std::ostream &out, Arg &&arg, Args &&... args)
{
  //  if constexpr (debugOutput) // probably will not get away with c++17
  if (debugOutput)
  {
    for (int i = 0; i < PCU_Comm_Peers(h); i++)
    {
      if (i == PCU_Comm_Self(h))
      {
        out << "Rank [" << i << "] " << std::forward<Arg>(arg);
        //((out << ", " << std::forward<Args>(args)), ...); // probably will not get away with c++17
        using expander = int[];
        (void)expander{0, (void(out << ", " << std::forward<Args>(args)), 0)...};
        out << "\n";
        out << std::flush;
      }
      // removed this: can't guarantee this is collectively called, order not ensured therefore.
      //PCU_Barrier();
    }
  }
}

template <typename... Args>
void Kill(PCU_t h, const int fid, Args &&... args)
{
  DebugParallelPrinter(h, std::cout, "***** CGNS ERROR", args...);

  if (PCU_Comm_Initialized(h))
  {
    cgp_close(fid);
    // Finalize the MPI environment.
    PCU_Comm_Free(h);
    MPI_Finalize();
    cgp_error_exit();
    exit(EXIT_FAILURE);
  }
  else
  {
    cg_close(fid);
    cg_error_exit();
    exit(EXIT_FAILURE);
  }
}

void Kill(PCU_t h, const int fid)
{
  DebugParallelPrinter(h, std::cout, "***** CGNS ERROR");

  if (PCU_Comm_Initialized(h))
  {
    cgp_close(fid);
    // Finalize the MPI environment.
    PCU_Comm_Free(h);
    MPI_Finalize();
    cgp_error_exit();
    exit(EXIT_FAILURE);
  }
  else
  {
    cg_close(fid);
    cg_error_exit();
    exit(EXIT_FAILURE);
  }
}

auto ReadCGNSCoords(PCU_t h, int cgid, int base, int zone, int ncoords, int nverts, const std::vector<cgsize_t> &, const apf::GlobalToVert &globalToVert)
{
  // Read min required as defined by consecutive range
  // make one based as ReadElements makes zero based
  const cgsize_t lowest = globalToVert.begin()->first + 1;
  const cgsize_t highest = globalToVert.rbegin()->first + 1;
  DebugParallelPrinter(h, std::cout, "From globalToVert ", lowest, highest, nverts);

  cgsize_t range_min[3];
  range_min[0] = lowest;
  range_min[1] = lowest;
  range_min[2] = lowest;
  cgsize_t range_max[3];
  range_max[0] = highest;
  range_max[1] = highest;
  range_max[2] = highest;

  std::vector<std::vector<double>> ordinates(3);
  const cgsize_t numToRead = range_max[0] - range_min[0] + 1; // one based
  ordinates[0].resize(numToRead, 0.0);
  ordinates[1].resize(numToRead, 0.0);
  ordinates[2].resize(numToRead, 0.0);

  std::vector<std::string> coord_names(3);
  coord_names[0].resize(CGIO_MAX_NAME_LENGTH + 1, ' ');
  coord_names[1].resize(CGIO_MAX_NAME_LENGTH + 1, ' ');
  coord_names[2].resize(CGIO_MAX_NAME_LENGTH + 1, ' ');

  for (int d = 0; d < ncoords; d++)
  {
    CGNS_ENUMT(DataType_t)
    datatype = CGNS_ENUMV(DataTypeNull);

    if (cg_coord_info(cgid, base, zone, d + 1, &datatype, &coord_names[d][0]))
    {
      std::cout << __LINE__ << " CGNS is dead " << std::endl;
      Kill(h, cgid);
    }
    const auto coord_name = std::string(coord_names[d].c_str());
    //boost::algorithm::trim(coord_name); // can't be bothered including boost

    auto &ordinate = ordinates[d];
    if (cgp_coord_read_data(cgid, base, zone, d + 1, &range_min[0], &range_max[0],
                            ordinate.data()))
    {
      std::cout << __LINE__ << " CGNS is dead " << std::endl;
      Kill(h, cgid);
    }
  }
  // to be clear, indices passed back are global, zero based
  // should only return what's inside globalToVert, but above I have to read more....
  std::map<int, std::array<double, 3>> points;
  int counter = lowest;
  for (int i = 0; i < numToRead; i++)
  {
    int zeroBased = counter - 1; // remove as per the addition above
    auto iter = globalToVert.find(zeroBased);
    if (iter != globalToVert.end())
    {
      points.insert(std::make_pair(zeroBased, std::array<double, 3>{{ordinates[0][i], ordinates[1][i], ordinates[2][i]}}));
    }
    counter++;
  }

  return points;
}

void SimpleElementPartition(PCU_t h, std::vector<cgsize_t> &numberToReadPerProc, std::vector<cgsize_t> &startingIndex, int el_start /* one based */, int el_end, int numElements)
{
  numberToReadPerProc.resize(PCU_Comm_Peers(h), 0);
  const cgsize_t blockSize = static_cast<cgsize_t>(
      std::floor(static_cast<double>(numElements) /
                 static_cast<double>(PCU_Comm_Peers(h))));

  DebugParallelPrinter(h, std::cout, "BlockSize", blockSize, "numElements", numElements, "comms Size", PCU_Comm_Peers(h));

  cgsize_t counter = 0;
  if (blockSize == 0 && numElements > 0)
    numberToReadPerProc[0] = numElements;
  else
  {
    bool keepGoing = true;
    while (keepGoing)
    {
      for (int p = 0; p < PCU_Comm_Peers(h) && keepGoing; p++)
      {
        if (counter + blockSize <= numElements)
        {
          counter += blockSize;
          numberToReadPerProc[p] += blockSize;
        }
        else
        {
          const auto &delta = numElements - counter;
          counter += delta;
          numberToReadPerProc[p] += delta;
          keepGoing = false;
        }
      }
    }
  }
  DebugParallelPrinter(h, std::cout, "Sanity check: Counter", counter, "numElements", numElements);

  DebugParallelPrinter(h, std::cout, "numberToReadPerProc for rank", PCU_Comm_Self(h), "is:", numberToReadPerProc[PCU_Comm_Self(h)]);

  startingIndex.resize(PCU_Comm_Peers(h), 0);
  startingIndex[0] = el_start;
  for (std::size_t i = 1; i < numberToReadPerProc.size(); i++)
  {
    if (numberToReadPerProc[i] > 0)
      startingIndex[i] = startingIndex[i - 1] + numberToReadPerProc[i - 1];
  }

  DebugParallelPrinter(h, std::cout, "Element start, end, numElements ", el_start,
                       el_end, numElements);

  DebugParallelPrinter(h, std::cout, "startingIndex for rank", PCU_Comm_Self(h), "is:", startingIndex[PCU_Comm_Self(h)]);

  DebugParallelPrinter(h, std::cout, "Returning from SimpleElementPartition \n");
}

using Pair = std::pair<cgsize_t, cgsize_t>;
using LocalElementRanges = std::vector<Pair>; // one based

auto ReadElements(PCU_t h, int cgid, int base, int zone, int section, int el_start /* one based */, int el_end, int numElements, int verticesPerElement, LocalElementRanges &localElementRanges)
{
  std::vector<cgsize_t> numberToReadPerProc;
  std::vector<cgsize_t> startingIndex;
  SimpleElementPartition(h, numberToReadPerProc, startingIndex, el_start, el_end, numElements);

  std::vector<cgsize_t> vertexIDs;
  if (numberToReadPerProc[PCU_Comm_Self(h)] > 0)
    vertexIDs.resize(numberToReadPerProc[PCU_Comm_Self(h)] * verticesPerElement,
                     -1234567);

  const auto start = startingIndex[PCU_Comm_Self(h)];
  const auto end = startingIndex[PCU_Comm_Self(h)] + numberToReadPerProc[PCU_Comm_Self(h)] - 1;
  DebugParallelPrinter(h, std::cout, "Range", start, "to", end, numberToReadPerProc[PCU_Comm_Self(h)]);
  //
  cgp_elements_read_data(cgid, base, zone, section, start,
                         end, vertexIDs.data());

  if (numberToReadPerProc[PCU_Comm_Self(h)] > 0)
  {
    localElementRanges.push_back(std::make_pair(start, end));
  }

  // remove CGNS one-based offset
  // to be clear these are the node ids defining the elements, not element ids
  for (auto &i : vertexIDs)
  {
    i = i - 1;
  }

  return std::make_tuple(vertexIDs, numberToReadPerProc[PCU_Comm_Self(h)]);
}

struct CGNSBCMeta
{
  std::string bocoName;
  CGNS_ENUMT(BCType_t)
  bocoType = CGNS_ENUMV(BCTypeNull);
  CGNS_ENUMT(PointSetType_t)
  ptsetType = CGNS_ENUMV(PointSetTypeNull);
  cgsize_t npnts = -1;
  std::vector<int> normalindices;
  cgsize_t normalListSize = -1;
  CGNS_ENUMT(DataType_t)
  normalDataType = CGNS_ENUMV(DataTypeNull);
  int ndataset = -1;
  CGNS_ENUMT(GridLocation_t)
  location = CGNS_ENUMV(GridLocationNull);
  std::string locationName;
  cgsize_t minElementId = -1;
  cgsize_t maxElementId = -1;
  std::vector<cgsize_t> bcElementIds;

  void Info() const
  {
    std::cout << "BC named: " << bocoName << ", located on: " << locationName << std::endl;
    std::cout << "\tHas " << npnts << " elements on the bc stored as a " << cg_PointSetTypeName(ptsetType) << std::endl;
    std::cout << "\tThe min and max Element Ids for this bcs are: " << minElementId << " " << maxElementId << std::endl;
    if (debugOutput)
    //if constexpr (debugOutput)// probably will not get away with c++17
    {
      std::cout << "\tThe element Ids that are tagged with this bcs are: ";
      for (const auto &i : bcElementIds)
        std::cout << i << " ";
      std::cout << std::endl;
    }
  }
};

struct BCInfo
{
  std::string bcName;                     // user provided
  std::string cgnsLocation;               // for debug
  std::unordered_set<cgsize_t> vertexIds; // zero-based global to relate to GlobalToVert
  apf::MeshTag *tag = nullptr;
  apf::Field *field = nullptr; // debug, outputting to vtk

  void Clean(apf::Mesh2 *m)
  {
    m->removeField(field);
    apf::destroyField(field);
    apf::removeTagFromDimension(m, tag, 0);
    m->destroyTag(tag);
  }

  void TagVertices(const int cgid, apf::Mesh2 *m, apf::GlobalToVert &globalToVert)
  {
    // tm = "tag marker" for that bcName, fm = "field marker" for that bcName
    tag = m->createIntTag(("debug_tm_" + bcName).c_str(), 1); // 1 is size of tag
    apf::MeshEntity *elem = nullptr;
    apf::MeshIterator *it = m->begin(0);
    int vals[1];
    vals[0] = 0;
    while ((elem = m->iterate(it)))
      m->setIntTag(elem, tag, vals);
    m->end(it);

    for (const auto &v : vertexIds)
    {
      auto iter = globalToVert.find(v);
      if (iter != globalToVert.end())
      {
        apf::MeshEntity *elem = iter->second;
        vals[0] = 1;
        m->setIntTag(elem, tag, vals);
      }
      else
      {
        Kill(m->getPCU()->GetCHandle(), cgid, "GlobalToVert lookup problem", v);
      }
    }

    if (debugOutput)
    //if constexpr (debugOutput) // probably will not get away with c++17
    {
      // for debug output, tags aren't written to vtk...
      apf::MeshEntity *elem = nullptr;
      apf::MeshIterator *it = m->begin(0);
      field = apf::createFieldOn(m, ("debug_fm_" + bcName).c_str(), apf::SCALAR);

      int vals[1];
      double dval;
      while ((elem = m->iterate(it)))
      {
        m->getIntTag(elem, tag, vals);
        dval = vals[0];
        apf::setScalar(field, elem, 0, dval);
      }
      m->end(it);
    }

    // Notes:
    // I do not exchange the tag values (even if that can be done).
    // I'm assuming in parallel all partitions that need the vertex that
    // falls within a given group will mark that vertex accordingly.
    // I assume therefore that vertices on a processor boundary, are marked
    // by all procs that share it.
    // TODO: generate test that proves this works. Done: for quads tested 4-way
  }

  /**
   * Meanings/concept follows: https://cgns.github.io/CGNS_docs_current/sids/misc.html section 12.4
    +---------------+--------------+--------------+--------------+-------------------------+
    |               | GridLocation | GridLocation | GridLocation |      GridLocation       |
    +---------------+--------------+--------------+--------------+-------------------------+
    | CellDimension | Vertex       | EdgeCenter   | *FaceCenter  | CellCenter              |
    | 1             | vertices     | -            | -            | cells (line elements)   |
    | 2             | vertices     | edges        | -            | cells (area elements)   |
    | 3             | vertices     | edges        | faces        | cells (volume elements) |
    +---------------+--------------+--------------+--------------+-------------------------+ 
  **/
  void TagBCEntities(const int cgid, apf::Mesh2 *m, apf::CGNSBCMap &cgnsBCMap)
  {
    const auto Add = [&cgnsBCMap](const std::string &location, const std::string &tagName, apf::MeshTag *field) {
      auto iter = cgnsBCMap.find(location);
      if (iter != cgnsBCMap.end())
      {
        apf::CGNSInfo info;
        info.cgnsBCSName = tagName;
        info.bcsMarkerTag = field;
        iter->second.push_back(info);
      }
      else
      {
        std::vector<apf::CGNSInfo> infos;
        apf::CGNSInfo info;
        info.cgnsBCSName = tagName;
        info.bcsMarkerTag = field;
        infos.push_back(info);
        cgnsBCMap.insert(std::make_pair(location, infos));
      }
    };

    const auto VertexLoop = [&m](const std::string &tagName, apf::MeshTag *vertexTag) {
      apf::MeshTag *bcTag = nullptr;
      bcTag = m->createIntTag(tagName.c_str(), 1); // 1 is size of tag

      apf::MeshIterator *vertIter = m->begin(0);
      apf::MeshEntity *vert = nullptr;
      int vals[1];
      vals[0] = 0;
      while ((vert = m->iterate(vertIter)))
      {
        m->setIntTag(vert, bcTag, vals);
      }
      m->end(vertIter);

      vertIter = m->begin(0);
      while ((vert = m->iterate(vertIter)))
      {
        bool allTagged = true;
        m->getIntTag(vert, vertexTag, vals);
        if (vals[0] == 0)
          allTagged = false;

        if (allTagged)
        {
          vals[0] = 1;
          m->setIntTag(vert, bcTag, vals);
        }
      }
      m->end(vertIter);

      if (debugOutput)
      //if constexpr (debugOutput) // probably will not get away with c++17
      {
        // for debug output, tags aren't written to vtk...
        apf::MeshEntity *elem = nullptr;
        apf::MeshIterator *it = m->begin(0);
        auto *field = apf::createFieldOn(m, ("debug_" + tagName).c_str(), apf::SCALAR);

        int vals[1];
        double dval;
        while ((elem = m->iterate(it)))
        {
          m->getIntTag(elem, bcTag, vals);
          dval = vals[0];
          apf::setScalar(field, elem, 0, dval);
        }
        m->end(it);
      }

      return bcTag;
    };

    const auto EdgeLoop = [&m](const std::string &tagName, apf::MeshTag *vertexTag) {
      apf::MeshTag *bcTag = nullptr;
      bcTag = m->createIntTag(tagName.c_str(), 1); // 1 is size of tag

      apf::MeshIterator *edgeIter = m->begin(1);
      apf::MeshEntity *edge = nullptr;
      int vals[1];
      vals[0] = 0;
      while ((edge = m->iterate(edgeIter)))
      {
        m->setIntTag(edge, bcTag, vals);
      }
      m->end(edgeIter);

      apf::Downward verts;
      edgeIter = m->begin(1);
      while ((edge = m->iterate(edgeIter)))
      {
        const auto numVerts = m->getDownward(edge, 0, verts);
        bool allTagged = true;
        for (int i = 0; i < numVerts; i++)
        {
          m->getIntTag(verts[i], vertexTag, vals);
          if (vals[0] == 0)
            allTagged = false;
        }
        if (allTagged)
        {
          vals[0] = 1;
          m->setIntTag(edge, bcTag, vals);
        }
      }
      m->end(edgeIter);

      if (debugOutput)
      //if constexpr (debugOutput) // probably will not get away with c++17
      {
        // for debug output, tags aren't written to vtk...
        apf::MeshEntity *elem = nullptr;
        apf::MeshIterator *it = m->begin(1);
        auto *field = apf::createField(m, ("debug_" + tagName).c_str(), apf::SCALAR, apf::getConstant(1));

        int vals[1];
        double dval;
        while ((elem = m->iterate(it)))
        {
          m->getIntTag(elem, bcTag, vals);
          dval = vals[0];
          apf::setScalar(field, elem, 0, dval);
        }
        m->end(it);
      }

      return bcTag;
    };

    const auto FaceLoop = [&m](const std::string &tagName, apf::MeshTag *vertexTag) {
      apf::MeshTag *bcTag = nullptr;
      bcTag = m->createIntTag(tagName.c_str(), 1); // 1 is size of tag

      apf::MeshIterator *faceIter = m->begin(2);
      apf::MeshEntity *face = nullptr;
      int vals[1];
      vals[0] = 0;
      while ((face = m->iterate(faceIter)))
      {
        m->setIntTag(face, bcTag, vals);
      }
      m->end(faceIter);

      apf::Downward verts;
      faceIter = m->begin(2);
      while ((face = m->iterate(faceIter)))
      {
        const auto numVerts = m->getDownward(face, 0, verts);
        bool allTagged = true;
        for (int i = 0; i < numVerts; i++)
        {
          m->getIntTag(verts[i], vertexTag, vals);
          if (vals[0] == 0)
            allTagged = false;
        }
        if (allTagged)
        {
          vals[0] = 1;
          m->setIntTag(face, bcTag, vals);
        }
      }
      m->end(faceIter);

      if (debugOutput)
      //if constexpr (debugOutput) // probably will not get away with c++17
      {
        // for debug output, tags aren't written to vtk...
        apf::MeshEntity *elem = nullptr;
        apf::MeshIterator *it = m->begin(2);
        auto *field = apf::createField(m, ("debug_" + tagName).c_str(), apf::SCALAR, apf::getConstant(2));

        int vals[1];
        double dval;
        while ((elem = m->iterate(it)))
        {
          m->getIntTag(elem, bcTag, vals);
          dval = vals[0];
          apf::setScalar(field, elem, 0, dval);
        }
        m->end(it);
      }

      return bcTag;
    };

    const auto CellLoop = [&m](const std::string &tagName, apf::MeshTag *vertexTag, int dim) {
      apf::MeshTag *bcTag = nullptr;
      bcTag = m->createIntTag(tagName.c_str(), 1); // 1 is size of tag

      apf::MeshIterator *cellIter = m->begin(dim);
      apf::MeshEntity *cell = nullptr;
      int vals[1];
      vals[0] = 0;
      while ((cell = m->iterate(cellIter)))
      {
        m->setIntTag(cell, bcTag, vals);
      }
      m->end(cellIter);

      apf::Downward verts;
      cellIter = m->begin(dim);
      while ((cell = m->iterate(cellIter)))
      {
        const auto numVerts = m->getDownward(cell, 0, verts);
        bool allTagged = true;
        for (int i = 0; i < numVerts; i++)
        {
          m->getIntTag(verts[i], vertexTag, vals);
          if (vals[0] == 0)
            allTagged = false;
        }
        if (allTagged)
        {
          vals[0] = 1;
          m->setIntTag(cell, bcTag, vals);
        }
      }
      m->end(cellIter);

      if (debugOutput)
      //if constexpr (debugOutput) // probably will not get away with c++17
      {
        // for debug output, tags aren't written to vtk...
        apf::MeshEntity *elem = nullptr;
        apf::MeshIterator *it = m->begin(dim);
        auto *field = apf::createField(m, ("debug_" + tagName).c_str(), apf::SCALAR, apf::getConstant(dim));

        int vals[1];
        double dval;
        while ((elem = m->iterate(it)))
        {
          m->getIntTag(elem, bcTag, vals);
          dval = vals[0];
          apf::setScalar(field, elem, 0, dval);
        }
        m->end(it);
      }

      return bcTag;
    };

    if (m->getDimension() == 3) // working with a 3D mesh
    {
      if (cgnsLocation == "Vertex")
      {
        apf::MeshTag *bcTag = VertexLoop(bcName, tag);
        Add("Vertex", bcName, bcTag);
      }
      else if (cgnsLocation == "EdgeCenter")
      {
        apf::MeshTag *bcTag = EdgeLoop(bcName, tag);
        Add("EdgeCenter", bcName, bcTag);
      }
      else if (cgnsLocation == "FaceCenter")
      {
        apf::MeshTag *bcTag = FaceLoop(bcName, tag);
        Add("FaceCenter", bcName, bcTag);
      }
      else if (cgnsLocation == "CellCenter")
      {
        apf::MeshTag *bcTag = CellLoop(bcName, tag, m->getDimension());
        Add("CellCenter", bcName, bcTag);
      }
      else
        Kill(m->getPCU()->GetCHandle(), cgid, "Unknown BC Type", cgnsLocation);
    }
    else if (m->getDimension() == 2) // working with a 2D mesh
    {
      if (cgnsLocation == "Vertex")
      {
        apf::MeshTag *bcTag = VertexLoop(bcName, tag);
        Add("Vertex", bcName, bcTag);
      }
      else if (cgnsLocation == "EdgeCenter")
      {
        apf::MeshTag *bcTag = EdgeLoop(bcName, tag);
        Add("EdgeCenter", bcName, bcTag);
      }
      else if (cgnsLocation == "FaceCenter")
      {
        PCU_ALWAYS_ASSERT_VERBOSE(true == false, "Can't have a FaceCenter BC in a 2D mesh");
      }
      else if (cgnsLocation == "CellCenter")
      {
        apf::MeshTag *bcTag = CellLoop(bcName, tag, m->getDimension());
        Add("CellCenter", bcName, bcTag);
      }
    }
    else if (m->getDimension() == 1) // working with a 1D mesh
    {
      if (cgnsLocation == "Vertex")
      {
        apf::MeshTag *bcTag = VertexLoop(bcName, tag);
        Add("Vertex", bcName, bcTag);
      }
      else if (cgnsLocation == "EdgeCenter")
      {
        PCU_ALWAYS_ASSERT_VERBOSE(true == false, "Can't have a EdgeCenter BC in a 1D mesh");
      }
      else if (cgnsLocation == "FaceCenter")
      {
        PCU_ALWAYS_ASSERT_VERBOSE(true == false, "Can't have a FaceCenter BC in a 1D mesh");
      }
      else if (cgnsLocation == "CellCenter")
      {
        apf::MeshTag *bcTag = CellLoop(bcName, tag, m->getDimension());
        Add("CellCenter", bcName, bcTag);
      }
    }

    if (!debugOutput)
      Clean(m);
  }
}; // namespace

void ReadBCInfo(PCU_t h, const int cgid, const int base, const int zone, const int nBocos, const int physDim, const int cellDim, const int nsections, std::vector<BCInfo> &bcInfos, const apf::GlobalToVert &globalToVert)
{
  // Read the BCS.
  std::vector<CGNSBCMeta> bcMetas(nBocos);
  bcInfos.resize(nBocos);
  //
  for (int boco = 1; boco <= nBocos; boco++)
  {
    auto &bcInfo = bcInfos[boco - 1];
    auto &bcMeta = bcMetas[boco - 1];
    bcMeta.normalindices.resize(physDim);

    bcMeta.bocoName.resize(CGIO_MAX_NAME_LENGTH + 1, ' ');
    bool pointRange = false;

    if (cg_boco_info(cgid, base, zone, boco, &bcMeta.bocoName[0], &bcMeta.bocoType,
                     &bcMeta.ptsetType, &bcMeta.npnts, bcMeta.normalindices.data(), &bcMeta.normalListSize,
                     &bcMeta.normalDataType, &bcMeta.ndataset))
      Kill(h, cgid, "Failed cg_boco_info");

    if (bcMeta.ptsetType == CGNS_ENUMV(PointList) || (bcMeta.ptsetType == CGNS_ENUMV(PointRange)))
    {
      bcMeta.bocoName = std::string(bcMeta.bocoName.c_str());
      //boost::algorithm::trim(bcMeta.bocoName); // can't be bothered including boost

      if (cg_boco_gridlocation_read(cgid, base, zone, boco, &bcMeta.location))
        Kill(h, cgid, "Failed cg_boco_gridlocation_read");

      bcMeta.locationName = cg_GridLocationName(bcMeta.location);

      if (bcMeta.ptsetType == CGNS_ENUMV(PointRange))
        pointRange = true;
    }
    else
    {
      Kill(h, cgid,
           "TODO: Can only work with "
           "PointList and PointRange BC "
           "types at the moment");
    }

    bcMeta.minElementId = -1;
    bcMeta.maxElementId = -1;
    bcMeta.bcElementIds.resize(bcMeta.npnts, -1);

    // here I read say elements with bcs as: 5, 3, 7, 9, and then read below ALL elements from 3->9,
    // but I only need, 3, 5, 7, 9, so don't need elements 4, 6, 8
    {
      if (bcMeta.locationName == "Vertex") // && cellDim == 1)
      {
        if (cg_boco_read(cgid, base, zone, boco, bcMeta.bcElementIds.data(), NULL))
          Kill(h, cgid, "Failed cg_boco_read");
      }
      else if (bcMeta.locationName == "EdgeCenter") // && cellDim == 2)
      {
        if (cg_boco_read(cgid, base, zone, boco, bcMeta.bcElementIds.data(), NULL))
          Kill(h, cgid, "Failed cg_boco_read");
      }
      else if (bcMeta.locationName == "FaceCenter") // && cellDim == 3)
      {
        if (cg_boco_read(cgid, base, zone, boco, bcMeta.bcElementIds.data(), NULL))
          Kill(h, cgid, "Failed cg_boco_read");
      }
      else if (bcMeta.locationName == "CellCenter") // && cellDim == 3)
      {
        if (cg_boco_read(cgid, base, zone, boco, bcMeta.bcElementIds.data(), NULL))
          Kill(h, cgid, "Failed cg_boco_read");
      }
      else
        Kill(h, cgid, "Failed Location test for BC Type", bcMeta.locationName,
             cellDim);

      if (pointRange)
      {
        // Check this is correct, I'm just trying to fill a contiguous range from [start, end]
        PCU_ALWAYS_ASSERT_VERBOSE(bcMeta.bcElementIds.size() == 2, "wrong size");
        const auto start = bcMeta.bcElementIds[0];
        const auto end = bcMeta.bcElementIds[1];
        const auto size = end - start + 1;
        bcMeta.bcElementIds.resize(size, -1);
        std::iota(std::begin(bcMeta.bcElementIds), std::end(bcMeta.bcElementIds), start);
      }

      bcMeta.minElementId = *std::min_element(bcMeta.bcElementIds.begin(), bcMeta.bcElementIds.end());
      bcMeta.maxElementId = *std::max_element(bcMeta.bcElementIds.begin(), bcMeta.bcElementIds.end());
      bcMeta.Info();
      bcInfo.bcName = bcMeta.bocoName;
      bcInfo.cgnsLocation = bcMeta.locationName;
    }

    std::vector<cgsize_t> vertexIDs;
    if (bcMeta.locationName != "Vertex")
    {
      std::unordered_set<cgsize_t> elementsToConsider;
      for (const auto &j : bcMeta.bcElementIds)
        elementsToConsider.insert(j);

      for (int section = 1; section <= nsections; section++)
      {
        std::string sectionName;
        sectionName.resize(CGIO_MAX_NAME_LENGTH + 1, ' ');

        CGNS_ENUMT(ElementType_t)
        elementType = CGNS_ENUMV(ElementTypeNull);
        cgsize_t el_start = -1;
        cgsize_t el_end = -1;
        int num_bndry = -1;
        int parent_flag = -1;
        //cgsize_t numElements = -1;
        int verticesPerElement = -1;

        cg_section_read(cgid, base, zone, section, &sectionName[0],
                        &elementType, &el_start, &el_end,
                        &num_bndry, &parent_flag);

        //numElements = el_end - el_start + 1;

        cg_npe(elementType, &verticesPerElement);

        const cgsize_t es = el_start;
        const cgsize_t ee = el_end;

        cgsize_t range_min = bcMeta.minElementId;
        cgsize_t range_max = bcMeta.maxElementId;

        bool doRead = true;
        if (es < bcMeta.minElementId && ee < bcMeta.minElementId)
          doRead = false;
        if (es > bcMeta.maxElementId && es > bcMeta.maxElementId)
          doRead = false;

        if (doRead)
        {
          if (es < bcMeta.minElementId && ee < bcMeta.minElementId) // clip lower
            range_max = ee;

          if (es < bcMeta.maxElementId && ee > bcMeta.maxElementId) // clip higher
            range_min = es;

          range_min = std::max(es, bcMeta.minElementId);
          range_max = std::min(ee, bcMeta.maxElementId);
          const cgsize_t range_num = range_max - range_min + 1;

          vertexIDs.resize(range_num * verticesPerElement,
                           -1234567);

          cg_elements_partial_read(cgid, base, zone, section, range_min,
                                   range_max, vertexIDs.data(), nullptr);

          cgsize_t counter = range_min;

          for (size_t i = 0; i < vertexIDs.size(); i += verticesPerElement)
          {
            if (elementsToConsider.count(counter))
            {
              auto last = std::min(vertexIDs.size(), i + verticesPerElement);
              std::vector<cgsize_t> vec(vertexIDs.begin() + i, vertexIDs.begin() + last);
              // std::cout << counter << " ";
              // for (const auto &j : vec)
              //   std::cout << j << " ";
              // std::cout << std::endl;

              for (const auto &k : vec)
              {
                const auto zb = k - 1; // make zero-based
                auto iter = globalToVert.find(zb);
                if (iter != globalToVert.end())
                  bcInfo.vertexIds.insert(zb);
              }
            }
            counter++;
          }
        }
      }
    }
    else
    {
      vertexIDs = bcMeta.bcElementIds;
      // std::cout << "Vertex BCS: ";
      // for (const auto &j : vertexIDs)
      //   std::cout << j << " ";
      // std::cout << std::endl;
      for (const auto &k : vertexIDs)
      {
        const auto zb = k - 1; // make zero-based
        auto iter = globalToVert.find(zb);
        if (iter != globalToVert.end())
          bcInfo.vertexIds.insert(zb);
      }
    }
  }
}

apf::Mesh2 *DoIt(PCU_t h, gmi_model *g, const std::string &fname, apf::CGNSBCMap &cgnsBCMap, const std::vector<std::pair<std::string, std::string>> &readMeshData)
{
  static_assert(std::is_same<cgsize_t, int>::value, "cgsize_t not compiled as int");

  int cgid = -1;
  auto comm = PCU_Get_Comm(h);
  cgp_mpi_comm(comm);
  cgp_pio_mode(CGNS_ENUMV(CGP_INDEPENDENT));
  cgp_open(fname.c_str(), CGNS_ENUMV(CG_MODE_READ), &cgid);

  int nbases = -1;
  cg_nbases(cgid, &nbases);
  if (nbases > 1)
  {
    std::cout << "CGNS file has " << nbases << " nbases, only used base=1" << std::endl;
  }

  std::string basename;
  basename.resize(CGIO_MAX_NAME_LENGTH + 1, ' ');

  int cellDim = -1;
  int physDim = -1;
  const int base = 1;
  cg_base_read(cgid, base, &basename[0], &cellDim, &physDim);
  const int readDim = cellDim;

  // Salome cgns is a bit on the odd side: cellDim, physDim, ncoords are not always consistent
  apf::Mesh2 *mesh = apf::makeEmptyMdsMesh(g, cellDim, false, static_cast<pcu::PCU*>(h.ptr));
  apf::GlobalToVert globalToVert;

  LocalElementRanges localElementRanges;
  using NewElements = std::vector<apf::MeshEntity *>;
  std::vector<NewElements> localElements;

  int nzones = -1;
  cg_nzones(cgid, base, &nzones);
  basename = std::string(basename.c_str());

  int nfam = -1;
  cg_nfamilies(cgid, base, &nfam);

  std::vector<BCInfo> bcInfos;

  for (int zone = 1; zone <= nzones; ++zone)
  {
    std::string zoneName;
    zoneName.resize(CGIO_MAX_NAME_LENGTH + 1, ' ');

    /* 
			3D unstructured: NVertex, NCell3D, NBoundVertex
			2D unstructured: NVertex, NCell2D, NBoundVertex
		*/
    std::array<cgsize_t, 3> sizes;
    cg_zone_read(cgid, base, zone, &zoneName[0], sizes.data());

    zoneName = std::string(zoneName.c_str());

    int ngrids = -1;
    cg_ngrids(cgid, base, zone, &ngrids);
    if (ngrids > 1)
    {
      std::cout << __LINE__ << " CGNS is dead " << std::endl;
      Kill(h, cgid);
    }
    int ncoords = -1;
    cg_ncoords(cgid, base, zone, &ncoords);

    CGNS_ENUMT(ZoneType_t)
    zonetype = CGNS_ENUMV(ZoneTypeNull);
    cg_zone_type(cgid, base, zone, &zonetype);

    int nsections = -1;
    if (zonetype == CGNS_ENUMV(Unstructured))
    {
      cg_nsections(cgid, base, zone, &nsections);
    }
    else
    {
      std::cout << __LINE__ << " CGNS is dead " << std::endl;
      Kill(h, cgid);
    }

    int nBocos = -1;
    cg_nbocos(cgid, base, zone, &nBocos);

    for (int section = 1; section <= nsections; section++)
    {
      std::string sectionName;
      sectionName.resize(CGIO_MAX_NAME_LENGTH + 1, ' ');

      CGNS_ENUMT(ElementType_t)
      elementType = CGNS_ENUMV(ElementTypeNull);
      cgsize_t el_start = -1;
      cgsize_t el_end = -1;
      int num_bndry = -1;
      int parent_flag = -1;
      cgsize_t numElements = -1;
      int verticesPerElement = -1;

      cg_section_read(cgid, base, zone, section, &sectionName[0],
                      &elementType, &el_start, &el_end,
                      &num_bndry, &parent_flag);

      numElements = el_end - el_start + 1;

      cg_npe(elementType, &verticesPerElement);

      const auto readElementsAndVerts = [&](const apf::Mesh2::Type &type) {
        const auto &ret = ReadElements(h, cgid, base, zone, section, el_start, el_end, numElements, verticesPerElement, localElementRanges);
        if (std::get<1>(ret) > 0)
        {
          const std::vector<cgsize_t> vertexIDs = std::get<0>(ret);
          std::vector<long> vertexIDs_l(vertexIDs.begin(), vertexIDs.end());
          localElements.emplace_back(apf::assemble(mesh, vertexIDs_l.data(), std::get<1>(ret), type, globalToVert)); // corresponding finalize below
          const auto nverts = sizes[0];
          const auto ordinates = ReadCGNSCoords(h, cgid, base, zone, ncoords, nverts, vertexIDs, globalToVert);

          for (const auto &p : globalToVert)
          {
            const auto pp = ordinates.at(p.first);
            apf::Vector3 point(pp[0], pp[1], pp[2]);
            mesh->setPoint(p.second, 0, point);

            // Don't know why I wrote it like this...
            // auto iter = globalToVert.find(p.first);
            // if (iter != globalToVert.end())
            // {
            //   mesh->setPoint(iter->second, 0, point);
            // }
            // else
            // {
            //   Kill(cgid, "GlobalToVert lookup problem");
            // }
          }
        }
      };

      if (elementType == CGNS_ENUMV(BAR_2))
      {
        if (readDim == 1)
          readElementsAndVerts(apf::Mesh2::EDGE);
      }
      else if (elementType == CGNS_ENUMV(QUAD_4))
      {
        if (readDim == 2)
          readElementsAndVerts(apf::Mesh2::QUAD);
      }
      else if (elementType == CGNS_ENUMV(TRI_3))
      {
        if (readDim == 2)
          readElementsAndVerts(apf::Mesh2::TRIANGLE);
      }
      else if (elementType == CGNS_ENUMV(TETRA_4))
      {
        if (readDim == 3)
          readElementsAndVerts(apf::Mesh2::TET);
      }
      else if (elementType == CGNS_ENUMV(PYRA_5))
      {
        if (readDim == 3)
          readElementsAndVerts(apf::Mesh2::PYRAMID);
      }
      else if (elementType == CGNS_ENUMV(HEXA_8))
      {
        if (readDim == 3)
          readElementsAndVerts(apf::Mesh2::HEX);
      }
      else
      {
        std::cout << __LINE__ << " CGNS is dead "
                  << " " << SupportedCGNSElementTypeToString(elementType) << std::endl;
        Kill(h, cgid);
      }
    }

    if (nBocos > 0)
    {
      std::cout << std::endl;
      std::cout << "Attempting to read BCS info "
                << " " << nBocos << std::endl;
      ReadBCInfo(h, cgid, base, zone, nBocos, physDim, cellDim, nsections, bcInfos, globalToVert);
      std::cout << std::endl;
    }

    int nsols = -1;
    if (cg_nsols(cgid, base, zone, &nsols))
      Kill(h, cgid, "1, ", nsols);

    std::vector<MeshData> meshData;
    if (nsols > 0)
    {
      for (int ns = 1; ns <= nsols; ns++)
      {
        CGNS_ENUMT(GridLocation_t)
        location;
        char sname[33];
        if (cg_sol_info(cgid, base, zone, ns,
                        sname, &location))
          Kill(h, cgid, "2");

        int nflds = -1;
        if (cg_nfields(cgid, base, zone, ns, &nflds))
          Kill(h, cgid, "3");

        if (nflds > 0)
        {
          for (int f = 1; f <= nflds; f++)
          {
            CGNS_ENUMT(DataType_t)
            datatype;
            char name[33];
            if (cg_field_info(cgid, base, zone, ns, f,
                              &datatype, name))
              Kill(h, cgid, "4");

            //std::cout << sname << " " << name << " " << f << " " << ns << " " << cg_DataTypeName(datatype) << " " << cg_GridLocationName(location) << std::endl;
            meshData.push_back(MeshData(ns, f, datatype, location, name));
          }
        }
      }
    }

    const auto index = [](const std::string &name) {
      for (int i = 0; i < 9; i++) // check for vectors and matrices
      {
        // pattern matching component writer in apfCGNS.cc, patterns must be kept in-sync with this file
        const std::string end("_[" + std::to_string(i) + "]");
        const std::size_t found = name.find(end);
        if (found != std::string::npos) // does name contain a [x] where x[0,9)
        {
          return i;
        }
      }
      return -1;
    };

    // process meshData and find vectors (and matrices)
    const auto findMatch = [&meshData, &index, &cgid](MeshData &other, std::vector<MeshDataGroup> &meshDataGroups) {
      MeshDataGroup group;
      for (auto &md : meshData)
      {
        if (md != other) // don't compare with yourself
        {
          if (md.process) // make sure md has not be marked as a component already
          {
            if (md.name.size() == other.name.size())
            {
              for (int i = 0; i < 9; i++) // check for vectors and matrices
              {
                // pattern matching component writer in apfCGNS.cc, patterns must be kept in-sync with this file
                const std::string end("_[" + std::to_string(i) + "]");
                const std::size_t found = md.name.find(end);
                if (found != std::string::npos) // does md contain a [x] where x[0,9)
                {
                  const auto stub = md.name.substr(0, found);         // get the part of the string not with a component index
                  const auto stubOther = other.name.substr(0, found); // get the part of the string not with a component index

                  if (stub == stubOther)
                  {
                    if (md.location == other.location) // are at the same location in the mesh
                    {
                      if (md.datatype == other.datatype) // have same data type
                      {
                        md.process = false;
                        other.process = false;
                        const auto otherIndex = index(other.name);
                        if (otherIndex != -1)
                          group.insert(otherIndex, other);
                        else
                          Kill(h, cgid, "Bad fail");
                        group.insert(i, md);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

      if (group.size() > 0)
        meshDataGroups.push_back(group);
    };

    // If a field or tag is already present, don't over-write or try to add
    // Can modify this later for different behaviour, but useful default.
    std::vector<std::string> existingNames;
    {
      apf::DynamicArray<apf::MeshTag *> tags;
      mesh->getTags(tags);
      for (std::size_t i = 0; i < tags.getSize(); ++i)
      {
        apf::MeshTag *t = tags[i];
        const std::string tagName(mesh->getTagName(t));
        existingNames.push_back(tagName);
      }
      for (int i = 0; i < mesh->countFields(); ++i)
      {
        apf::Field *f = mesh->getField(i);
        const std::string fieldName(f->getName());
        existingNames.push_back(fieldName);
      }
    }

    for (auto &md : meshData)
    {
      if (md.process)
      {
        for (const auto &n : existingNames)
        {
          if (n == md.name)
            md.process = false;
        }
      }
    }

    std::vector<MeshDataGroup> meshDataGroups;
    for (auto &md : meshData)
    {
      if (md.process)
        findMatch(md, meshDataGroups);
    }

    for (auto &md : meshData)
    {
      if (md.process)
      {
        MeshDataGroup group;
        group.insert(0, md);
        meshDataGroups.push_back(group);
        md.process = false;
      }
    }

    for (auto &md : meshDataGroups)
    {
      if (md.process)
      {
        bool read = false;
        for (const auto &n : readMeshData)
        {
          if (n.second == md.name())
          {
            if (n.first == cg_GridLocationName(md.location()))
            {
              std::cout << md.name() << " " << n.first << " " << n.second << " " << cg_GridLocationName(md.location()) << std::endl;
              read = true;
            }
          }
        }
        if (!read)
          md.process = false;
      }
    }

    // add all double's as fields, even though they once may have been tags, since
    // that info is lost in the cgns-writer (cos of cgns...)
    for (const auto &md : meshDataGroups)
    {
      if (md.process)
      {
        if (md.location() == CGNS_ENUMV(Vertex))
        {
          if (md.datatype() == CGNS_ENUMV(Integer))
          {
          }
          else if (md.datatype() == CGNS_ENUMV(RealDouble))
          {
            const cgsize_t lowest = globalToVert.begin()->first + 1;   // one based
            const cgsize_t highest = globalToVert.rbegin()->first + 1; // one based

            cgsize_t range_min[3];
            range_min[0] = lowest;
            range_min[1] = lowest;
            range_min[2] = lowest;
            cgsize_t range_max[3];
            range_max[0] = highest;
            range_max[1] = highest;
            range_max[2] = highest;
            const cgsize_t numToRead = range_max[0] - range_min[0] + 1; // one based

            apf::Field *field = nullptr;
            if (md.size() == 1)
              field = apf::createFieldOn(mesh, md.name().c_str(), apf::SCALAR);
            else if (md.size() == 3)
              field = apf::createFieldOn(mesh, md.name().c_str(), apf::VECTOR);
            else if (md.size() == 9)
              field = apf::createFieldOn(mesh, md.name().c_str(), apf::MATRIX);
            else
              Kill(h, cgid, "Tensor size not accounted for");

            double scalar(-123456.0);
            apf::Vector3 vector3(-123456.0, -123456.0, -123456.0);
            apf::Matrix3x3 matrix3x3(-123456.0, -123456.0, -123456.0, -123456.0, -123456.0, -123456.0, -123456.0, -123456.0, -123456.0);

            apf::MeshEntity *elem = nullptr;
            apf::MeshIterator *it = mesh->begin(0);
            while ((elem = mesh->iterate(it)))
            {
              if (md.size() == 1)
                apf::setScalar(field, elem, 0, scalar);
              else if (md.size() == 3)
                apf::setVector(field, elem, 0, vector3);
              else if (md.size() == 9)
                apf::setMatrix(field, elem, 0, matrix3x3);
              else
                Kill(h, cgid, "Tensor size not accounted for");
            }
            mesh->end(it);

            using CGNSType = double;
            std::vector<CGNSType> meshVals(numToRead, -123456.0);
            for (std::size_t i = 0; i < md.size(); i++)
            {
              if (cgp_field_read_data(cgid, base, zone, md.sIndex(i), md.fIndex(i),
                                      &range_min[0], &range_max[0], meshVals.data()))
                Kill(h, cgid, "Failed cgp_field_read_data");

              cgsize_t counter = lowest;
              for (cgsize_t it = 0; it < numToRead; it++)
              {
                cgsize_t zeroBased = counter - 1; // remove as per the addition above
                auto iter = globalToVert.find(zeroBased);
                if (iter != globalToVert.end())
                {
                  auto *elem = iter->second;
                  if (md.size() == 1)
                  {
                    apf::setScalar(field, elem, 0, meshVals.at(it));
                  }
                  else if (md.size() == 3)
                  {
                    apf::getVector(field, elem, 0, vector3);
                    vector3[i] = meshVals.at(it);
                    apf::setVector(field, elem, 0, vector3);
                  }
                  else if (md.size() == 9)
                  {
                    apf::getMatrix(field, elem, 0, matrix3x3);
                    if (i < 3)
                      matrix3x3[0][i] = meshVals.at(it);
                    else if (i > 2 && i < 6)
                      matrix3x3[1][i - 3] = meshVals.at(it);
                    else if (i > 5 && i < 9)
                      matrix3x3[2][i - 6] = meshVals.at(it);
                    apf::setMatrix(field, elem, 0, matrix3x3);
                  }
                  else
                    Kill(h, cgid, "Tensor size not accounted for");
                }
                counter++;
              }
            }
          }
          else if (md.datatype() == CGNS_ENUMV(LongInteger))
          {
          }
          else
          {
            Kill(h, cgid, "Don't know how to process this at the moment");
          }
        }
        else if (md.location() == CGNS_ENUMV(CellCenter))
        {
          if (md.datatype() == CGNS_ENUMV(Integer))
          {
          }
          else if (md.datatype() == CGNS_ENUMV(RealDouble))
          {
            // HUGE assumption here, I assume the order the cgns elements are read (and therefore their global number in the file)
            // is the same order they are created and added to the mesh by buildElements.
            // I combine localElementRanges to give me the global indices, with vector<NewElements> and hope it works
            PCU_ALWAYS_ASSERT_VERBOSE(localElementRanges.size() == localElements.size(),
                                      "Size don't match for element/number ranges");

            const auto dim = mesh->getDimension();
            apf::Field *field = nullptr;
            md.info();
            if (md.size() == 1)
              field = apf::createField(mesh, md.name().c_str(), apf::SCALAR, apf::getConstant(dim));
            else if (md.size() == 3)
              field = apf::createField(mesh, md.name().c_str(), apf::VECTOR, apf::getConstant(dim));
            else if (md.size() == 9)
              field = apf::createField(mesh, md.name().c_str(), apf::MATRIX, apf::getConstant(dim));
            else
              Kill(h, cgid, "Tensor size not accounted for");

            double scalar(-123456.0);
            apf::Vector3 vector3(-123456.0, -123456.0, -123456.0);
            apf::Matrix3x3 matrix3x3(-123456.0, -123456.0, -123456.0, -123456.0, -123456.0, -123456.0, -123456.0, -123456.0, -123456.0);

            apf::MeshEntity *elem = nullptr;
            apf::MeshIterator *it = mesh->begin(dim);
            while ((elem = mesh->iterate(it)))
            {
              if (md.size() == 1)
                apf::setScalar(field, elem, 0, scalar);
              else if (md.size() == 3)
                apf::setVector(field, elem, 0, vector3);
              else if (md.size() == 9)
                apf::setMatrix(field, elem, 0, matrix3x3);
              else
                Kill(h, cgid, "Tensor size not accounted for");
            }
            mesh->end(it);

            for (std::size_t r = 0; r < localElementRanges.size(); r++)
            {
              const auto range = localElementRanges[r];
              const cgsize_t lowest = range.first;   // one based
              const cgsize_t highest = range.second; // one based

              cgsize_t range_min[3];
              range_min[0] = lowest;
              range_min[1] = lowest;
              range_min[2] = lowest;
              cgsize_t range_max[3];
              range_max[0] = highest;
              range_max[1] = highest;
              range_max[2] = highest;
              const std::size_t numToRead = range_max[0] - range_min[0] + 1; // one based
              PCU_ALWAYS_ASSERT_VERBOSE(numToRead == localElements[r].size(),
                                        "Size don't match for element/number sub-ranges");

              using CGNSType = double;
              std::vector<CGNSType> meshVals(numToRead, -123456.0);
              for (std::size_t i = 0; i < md.size(); i++)
              {
                if (cgp_field_read_data(cgid, base, zone, md.sIndex(i), md.fIndex(i),
                                        &range_min[0], &range_max[0], meshVals.data()))
                  Kill(h, cgid, "Failed cgp_field_read_data");

                for (std::size_t it = 0; it < numToRead; it++)
                {
                  elem = localElements[r][it];

                  if (md.size() == 1)
                  {
                    apf::setScalar(field, elem, 0, meshVals.at(it));
                  }
                  else if (md.size() == 3)
                  {
                    apf::getVector(field, elem, 0, vector3);
                    vector3[i] = meshVals.at(it);
                    apf::setVector(field, elem, 0, vector3);
                  }
                  else if (md.size() == 9)
                  {
                    apf::getMatrix(field, elem, 0, matrix3x3);
                    if (i < 3)
                      matrix3x3[0][i] = meshVals.at(it);
                    else if (i > 2 && i < 6)
                      matrix3x3[1][i - 3] = meshVals.at(it);
                    else if (i > 5 && i < 9)
                      matrix3x3[2][i - 6] = meshVals.at(it);
                    apf::setMatrix(field, elem, 0, matrix3x3);
                  }
                  else
                  {
                    Kill(h, cgid, "Tensor size not accounted for");
                  }
                }
              }
            }
          }
          else if (md.datatype() == CGNS_ENUMV(LongInteger))
          {
          }
          else
          {
            Kill(h, cgid, "Don't know how to process this at the moment");
          }
        }
        else
        {
          Kill(h, cgid, "Don't know how to process this at the moment");
        }
      }
    }
  }

  // free up memory
  if (PCU_Comm_Initialized(h))
    cgp_close(cgid);
  else
    cg_close(cgid);

  {
    apf::MeshTag *tag = nullptr;
    tag = mesh->createIntTag("origCGNSGlobalVertID", 1);       // 1 is size of tag
    const cgsize_t lowest = globalToVert.begin()->first + 1;   // one based
    const cgsize_t highest = globalToVert.rbegin()->first + 1; // one based
    const cgsize_t numToRead = highest - lowest + 1;           // one based
    cgsize_t counter = lowest;
    for (cgsize_t it = 0; it < numToRead; it++)
    {
      cgsize_t zeroBased = counter - 1; // remove as per the addition above
      auto iter = globalToVert.find(zeroBased);
      if (iter != globalToVert.end())
      {
        const int oneBased = zeroBased + 1;
        mesh->setIntTag(iter->second, tag, &oneBased);
      }
      counter++;
    }
  }

  {
    apf::MeshTag *tag = nullptr;
    tag = mesh->createIntTag("origCGNSGlobalElemID", 1); // 1 is size of tag
    for (std::size_t r = 0; r < localElementRanges.size(); r++)
    {
      const auto range = localElementRanges[r];
      const cgsize_t lowest = range.first;                // one based
      const cgsize_t highest = range.second;              // one based
      const std::size_t numToRead = highest - lowest + 1; // one based
      int counter = lowest;
      for (std::size_t it = 0; it < numToRead; it++)
      {
        apf::MeshEntity *elem = nullptr;
        elem = localElements[r][it];
        mesh->setIntTag(elem, tag, &counter);
        counter++;
      }
    }
  }
  // not sure of the order of these three
  // and with reference to:

  apf::finalise(mesh, globalToVert);
  mesh->acceptChanges();
  apf::alignMdsRemotes(mesh);

  {
    apf::GlobalNumbering *gn = nullptr;
    gn = apf::makeGlobal(apf::numberOwnedNodes(mesh, "vert Idx"));
    apf::synchronize(gn);
  }
  // no synchronize call
  // https://github.com/SNLComputation/Albany/blob/master/src/disc/pumi/Albany_APFDiscretization.cpp @ various place throughout file
  // https://github.com/SCOREC/core/issues/249
  {
    apf::makeGlobal(apf::numberElements(mesh, "elem Idx"));
  }

  for (auto &bc : bcInfos)
  {
    bc.TagVertices(cgid, mesh, globalToVert);
  }

  for (auto &bc : bcInfos)
  {
    bc.TagBCEntities(cgid, mesh, cgnsBCMap);
  }

  apf::deriveMdsModel(mesh);
  apf::verify(mesh, true);

  return mesh;
} // namespace

apf::Mesh2 *DoIt(PCU_t h, gmi_model *g, const std::string &fname, apf::CGNSBCMap &cgnsBCMap)
{
  std::vector<std::pair<std::string, std::string>> meshData;
  return DoIt(h, g, fname, cgnsBCMap, meshData);
}

} // namespace

namespace apf
{

// caller needs to bring up and pull down mpi/pcu: mpi/pcu is required and assumed.
Mesh2 *loadMdsFromCGNS(gmi_model *g, const char *fname, apf::CGNSBCMap &cgnsBCMap, const std::vector<std::pair<std::string, std::string>> &meshData)
{
  PCU_t h = PCU_Get_Global_Handle();
  return loadMdsFromCGNS(h, g, fname, cgnsBCMap, meshData);
}

Mesh2 *loadMdsFromCGNS(PCU_t h, gmi_model *g, const char *fname, apf::CGNSBCMap &cgnsBCMap, const std::vector<std::pair<std::string, std::string>> &meshData)
{
#ifdef HAVE_CGNS
  Mesh2 *m = DoIt(h, g, fname, cgnsBCMap, meshData);
  return m;
#else
  Mesh2 *m = nullptr;
  PCU_ALWAYS_ASSERT_VERBOSE(m != nullptr,
                            "Build with ENABLE_CGNS to allow this functionality.");
  exit(EXIT_FAILURE);
  return m;
#endif
}

// caller needs to bring up and pull down mpi/pcu: mpi/pcu is required and assumed.
Mesh2 *loadMdsFromCGNS(gmi_model *g, const char *fname, apf::CGNSBCMap &cgnsBCMap)
{
  PCU_t h = PCU_Get_Global_Handle();
  return loadMdsFromCGNS(h, g, fname, cgnsBCMap);
}

Mesh2 *loadMdsFromCGNS(PCU_t h, gmi_model *g, const char *fname, apf::CGNSBCMap &cgnsBCMap)
{
#ifdef HAVE_CGNS
  Mesh2 *m = DoIt(h, g, fname, cgnsBCMap);
  return m;
#else
  Mesh2 *m = nullptr;
  PCU_ALWAYS_ASSERT_VERBOSE(m != nullptr,
                            "Build with ENABLE_CGNS to allow this functionality.");
  exit(EXIT_FAILURE);
  return m;
#endif
}

} // namespace apf
