/* Author: Dr Andrew Parker (2019) - FGE Ltd 

   Notes: these functions have been exclusively developed using meshes from Salome

*/

#include "apf.h"
#include "apfMDS.h"
#include "apfMesh2.h"
#include "apfShape.h"
#include "gmi.h" /* this is for gmi_getline... */
#include <lionPrint.h>
#include <PCU.h>
#include <apf.h>
#include <apfConvert.h>
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
#include <vector>
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
  else if (elementType == CGNS_ENUMV(PENTA_6))
    return "PENTA_6";
  else
  {
    std::cout << __LINE__ << " "
              << "No known element type "
              << " " << cg_ElementTypeName(elementType) << std::endl;
    return "";
  }
}

template <typename Arg, typename... Args>
void DebugParallelPrinter(std::ostream &out, Arg &&arg, Args &&... args)
{
  //  if constexpr (debugOutput) // probably will not get away with c++17
  if (debugOutput)
  {
    for (int i = 0; i < PCU_Comm_Peers(); i++)
    {
      if (i == PCU_Comm_Self())
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
void Kill(const int fid, Args &&... args)
{
  DebugParallelPrinter(std::cout, "***** CGNS ERROR", args...);

  if (PCU_Comm_Initialized())
  {
    cgp_close(fid);
    // Finalize the MPI environment.
    PCU_Comm_Free();
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

void Kill(const int fid)
{
  if (PCU_Comm_Initialized())
  {
    cgp_close(fid);
    // Finalize the MPI environment.
    PCU_Comm_Free();
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

auto ReadCGNSCoords(int cgid, int base, int zone, int ncoords, int nverts, const std::vector<cgsize_t> &, const apf::GlobalToVert &globalToVert)
{
  // Read min required as defined by consecutive range
  // make one based as ReadElements makes zero based
  const int lowest = globalToVert.begin()->first + 1;
  const int highest = globalToVert.rbegin()->first + 1;
  DebugParallelPrinter(std::cout, "From globalToVert ", lowest, highest, nverts);

  cgsize_t range_min[3];
  range_min[0] = lowest;
  range_min[1] = lowest;
  range_min[2] = lowest;
  cgsize_t range_max[3];
  range_max[0] = highest;
  range_max[1] = highest;
  range_max[2] = highest;

  std::vector<std::vector<double>> ordinates(3);
  const auto numToRead = range_max[0] - range_min[0] + 1;
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
      Kill(cgid);
    }
    const auto coord_name = std::string(coord_names[d].c_str());
    //boost::algorithm::trim(coord_name); // can't be bothered including boost

    auto &ordinate = ordinates[d];
    if (cgp_coord_read_data(cgid, base, zone, d + 1, &range_min[0], &range_max[0],
                            ordinate.data()))
    {
      std::cout << __LINE__ << " CGNS is dead " << std::endl;
      Kill(cgid);
    }
  }
  // to be clear, indices passed back are global, zero based
  std::map<int, std::array<double, 3>> points;
  int counter = lowest;
  for (int i = 0; i < numToRead; i++)
  {
    int zeroBased = counter - 1;
    points.insert(std::make_pair(zeroBased, std::array<double, 3>{{ordinates[0][i], ordinates[1][i], ordinates[2][i]}}));
    counter++;
  }

  return points;
}

void SimpleElementPartition(std::vector<cgsize_t> &numberToReadPerProc, std::vector<cgsize_t> &startingIndex, int el_start /* one based */, int el_end, int numElements)
{
  numberToReadPerProc.resize(PCU_Comm_Peers(), 0);
  const cgsize_t blockSize = static_cast<cgsize_t>(
      std::floor(static_cast<double>(numElements) /
                 static_cast<double>(PCU_Comm_Peers())));

  DebugParallelPrinter(std::cout, "BlockSize", blockSize, "numElements", numElements, "comms Size", PCU_Comm_Peers());

  cgsize_t counter = 0;
  if (blockSize == 0 && numElements > 0)
    numberToReadPerProc[0] = numElements;
  else
  {
    bool keepGoing = true;
    while (keepGoing)
    {
      for (int p = 0; p < PCU_Comm_Peers() && keepGoing; p++)
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
  DebugParallelPrinter(std::cout, "Sanity check: Counter", counter, "numElements", numElements);

  DebugParallelPrinter(std::cout, "numberToReadPerProc for rank", PCU_Comm_Self(), "is:", numberToReadPerProc[PCU_Comm_Self()]);

  startingIndex.resize(PCU_Comm_Peers(), 0);
  startingIndex[0] = el_start;
  for (std::size_t i = 1; i < numberToReadPerProc.size(); i++)
  {
    if (numberToReadPerProc[i] > 0)
      startingIndex[i] = startingIndex[i - 1] + numberToReadPerProc[i - 1];
  }

  DebugParallelPrinter(std::cout, "Element start, end, numElements ", el_start,
                       el_end, numElements);

  DebugParallelPrinter(std::cout, "startingIndex for rank", PCU_Comm_Self(), "is:", startingIndex[PCU_Comm_Self()]);

  DebugParallelPrinter(std::cout, "Returning from SimpleElementPartition \n");
}

auto ReadElements(int cgid, int base, int zone, int section, int el_start /* one based */, int el_end, int numElements, int verticesPerElement)
{
  std::vector<cgsize_t> numberToReadPerProc;
  std::vector<cgsize_t> startingIndex;
  SimpleElementPartition(numberToReadPerProc, startingIndex, el_start, el_end, numElements);

  std::vector<cgsize_t> vertexIDs;
  if (numberToReadPerProc[PCU_Comm_Self()] > 0)
    vertexIDs.resize(numberToReadPerProc[PCU_Comm_Self()] * verticesPerElement,
                     -1234567);

  const auto start = startingIndex[PCU_Comm_Self()];
  const auto end = startingIndex[PCU_Comm_Self()] + numberToReadPerProc[PCU_Comm_Self()] - 1;
  DebugParallelPrinter(std::cout, "Range", start, "to", end, numberToReadPerProc[PCU_Comm_Self()]);
  //
  cgp_elements_read_data(cgid, base, zone, section, start,
                         end, vertexIDs.data());

  // remove CGNS one-based offset
  // to be clear these are the node ids defining the elements, not element ids
  for (auto &i : vertexIDs)
  {
    i = i - 1;
  }

  return std::make_tuple(vertexIDs, numberToReadPerProc[PCU_Comm_Self()]);
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
  apf::Field *field = nullptr;

  void Clean(apf::Mesh *m)
  {
    m->removeField(field);
    apf::destroyField(field);
    apf::removeTagFromDimension(m, tag, 0);
    m->destroyTag(tag);
  }

  void TagVertices(const int cgid, apf::Mesh *m, apf::GlobalToVert &globalToVert)
  {
    tag = m->createIntTag(bcName.c_str(), 1); // 1 is size of tag
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
        Kill(cgid, "GlobalToVert lookup problem", v);
      }
    }

    if (debugOutput)
    //if constexpr (debugOutput) // probably will not get away with c++17
    {
      // for debug output, tags aren't written to vtk...
      apf::MeshEntity *elem = nullptr;
      apf::MeshIterator *it = m->begin(0);
      field = apf::createFieldOn(m, bcName.c_str(), apf::SCALAR);

      int vals[1];
      double dval[1];
      while ((elem = m->iterate(it)))
      {
        m->getIntTag(elem, tag, vals);
        dval[0] = vals[0];
        apf::setComponents(field, elem, 0, dval);
      }
      m->end(it);
    }

    // Notes:
    // I do not exchange the tag values (even if that can be done).
    // I'm assuming in parallel all partition that need the vertex that
    // falls within a given group mark that vertex accordingly.
    // I assume therefore that vertices on a processor boundary, are marked
    // by all procs that share it.
    // TODO: generate test that proves this works
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
  void TagBCEntities(const int cgid, apf::Mesh *m, apf::CGNSBCMap &cgnsBCMap)
  {
    if (m->getDimension() == 3) // working with a 3D mesh
    {
      if (cgnsLocation == "Vertex")
      {
        apf::Field *field = nullptr;
        const std::string fieldName = "VertBC_" + bcName;
        field = apf::createFieldOn(m, fieldName.c_str(), apf::SCALAR);

        double dval[1];
        dval[0] = 0.0;
        apf::MeshIterator *vertIter = m->begin(0);
        apf::MeshEntity *vert = nullptr;
        while ((vert = m->iterate(vertIter)))
        {
          apf::setComponents(field, vert, 0, dval);
        }
        m->end(vertIter);

        int vals[1];
        vertIter = m->begin(0);
        while ((vert = m->iterate(vertIter)))
        {
          bool allTagged = true;
          m->getIntTag(vert, tag, vals);
          if (vals[0] == 0)
            allTagged = false;

          if (allTagged)
          {
            std::cout << "Flagged vertices " << 1 << " " << allTagged << std::endl;
            dval[0] = 1.0;
            apf::setComponents(field, vert, 0, dval);
          }
        }
        auto iter = cgnsBCMap["Vertex"];
        iter.push_back(std::make_pair(fieldName, field));
      }
      else if (cgnsLocation == "EdgeCenter")
      {
        apf::Field *field = nullptr;
        const std::string fieldName = "EdgeBC_" + bcName;
        field = apf::createField(m, fieldName.c_str(), apf::SCALAR, apf::getConstant(1));

        double dval[1];
        dval[0] = 0.0;
        apf::MeshIterator *edgeIter = m->begin(1);
        apf::MeshEntity *edge = nullptr;
        while ((edge = m->iterate(edgeIter)))
        {
          apf::setComponents(field, edge, 0, dval);
        }
        m->end(edgeIter);

        apf::Downward verts;
        int vals[1];
        edgeIter = m->begin(1);
        while ((edge = m->iterate(edgeIter)))
        {
          const auto numVerts = m->getDownward(edge, 0, verts);
          bool allTagged = true;
          for (int i = 0; i < numVerts; i++)
          {
            m->getIntTag(verts[i], tag, vals);
            if (vals[0] == 0)
              allTagged = false;
          }
          if (allTagged)
          {
            //std::cout << "Flagged edges " << numVerts << " " << allTagged << std::endl;
            dval[0] = 1.0;
            apf::setComponents(field, edge, 0, dval);
          }
        }
        auto iter = cgnsBCMap["EdgeCenter"];
        iter.push_back(std::make_pair(fieldName, field));
      }
      else if (cgnsLocation == "FaceCenter")
      {
        apf::Field *field = nullptr;
        const std::string fieldName = "FaceBC_" + bcName;
        field = apf::createField(m, fieldName.c_str(), apf::SCALAR, apf::getConstant(2));

        double dval[1];
        dval[0] = 0.0;
        apf::MeshIterator *faceIter = m->begin(2);
        apf::MeshEntity *face = nullptr;
        while ((face = m->iterate(faceIter)))
        {
          apf::setComponents(field, face, 0, dval);
        }
        m->end(faceIter);

        apf::Downward verts;
        faceIter = m->begin(2);
        int vals[1];
        while ((face = m->iterate(faceIter)))
        {
          const auto numVerts = m->getDownward(face, 0, verts);
          bool allTagged = true;
          for (int i = 0; i < numVerts; i++)
          {
            m->getIntTag(verts[i], tag, vals);
            if (vals[0] == 0)
              allTagged = false;
          }
          if (allTagged)
          {
            //std::cout << "Flagged faces " << numVerts << " " << allTagged << std::endl;
            dval[0] = 1.0;
            apf::setComponents(field, face, 0, dval);
          }
        }
        auto iter = cgnsBCMap["FaceCenter"];
        iter.push_back(std::make_pair(fieldName, field));
      }
      else if (cgnsLocation == "CellCenter")
      {
        {
          apf::Field *field = nullptr;
          const std::string fieldName = "CellBC_" + bcName;
          field = apf::createField(m, fieldName.c_str(), apf::SCALAR, apf::getConstant(3));

          double dval[1];
          dval[0] = 0.0;
          apf::MeshIterator *cellIter = m->begin(3);
          apf::MeshEntity *cell = nullptr;
          while ((cell = m->iterate(cellIter)))
          {
            apf::setComponents(field, cell, 0, dval);
          }
          m->end(cellIter);

          apf::Downward verts;
          cellIter = m->begin(3);
          int vals[1];
          while ((cell = m->iterate(cellIter)))
          {
            const auto numVerts = m->getDownward(cell, 0, verts);
            bool allTagged = true;
            for (int i = 0; i < numVerts; i++)
            {
              m->getIntTag(verts[i], tag, vals);
              if (vals[0] == 0)
                allTagged = false;
            }
            if (allTagged)
            {
              //std::cout << "Flagged cells " << numVerts << " " << allTagged << std::endl;
              dval[0] = 1.0;
              apf::setComponents(field, cell, 0, dval);
            }
          }
          auto iter = cgnsBCMap["CellCenter"];
          iter.push_back(std::make_pair(fieldName, field));
        }
        // { // more verbose example of iterating the mesh
        //   apf::Field *field = nullptr;
        //   const std::string fieldName = "CellBC_other_" + bcName;
        //   field = apf::createField(m, fieldName.c_str(), apf::SCALAR, apf::getConstant(3));

        //   double dval[1];
        //   dval[0] = 0.0;
        //   apf::MeshIterator *cellIter = m->begin(3);
        //   apf::MeshEntity *cell = nullptr;
        //   while ((cell = m->iterate(cellIter)))
        //   {
        //     apf::setComponents(field, cell, 0, dval);
        //   }
        //   m->end(cellIter);

        //   apf::Downward faces;
        //   apf::Downward edges;
        //   apf::Downward verts;
        //   cellIter = m->begin(3);
        //   int vals[1];
        //   while ((cell = m->iterate(cellIter)))
        //   {
        //     bool allTagged = true;
        //     const auto numFaces = m->getDownward(cell, 2, faces);
        //     for (int f = 0; f < numFaces; f++)
        //     {
        //       const auto numEdges = m->getDownward(faces[f], 1, edges);
        //       for (int e = 0; e < numEdges; e++)
        //       {
        //         const auto numVerts = m->getDownward(edges[e], 0, verts);
        //         for (int i = 0; i < numVerts; i++)
        //         {
        //           m->getIntTag(verts[i], tag, vals);
        //           if (vals[0] == 0)
        //             allTagged = false;
        //         }
        //       }
        //     }
        //     if (allTagged)
        //     {
        //       std::cout << "Flagged cells " << allTagged << std::endl;
        //       dval[0] = 1.0;
        //       apf::setComponents(field, cell, 0, dval);
        //     }
        //   }
        //   auto iter = cgnsBCMap["CellCenter"];
        //   iter.push_back(std::make_pair(fieldName, field));
        // }
      }
      else
        Kill(cgid, "Unknown BC Type", cgnsLocation);
    }
    else if (m->getDimension() == 2) // working with a 2D mesh
    {
      if (cgnsLocation == "Vertex")
      {
        apf::Field *field = nullptr;
        const std::string fieldName = "VertBC_" + bcName;
        field = apf::createFieldOn(m, fieldName.c_str(), apf::SCALAR);

        double dval[1];
        dval[0] = 0.0;
        apf::MeshIterator *vertIter = m->begin(0);
        apf::MeshEntity *vert = nullptr;
        while ((vert = m->iterate(vertIter)))
        {
          apf::setComponents(field, vert, 0, dval);
        }
        m->end(vertIter);

        int vals[1];
        vertIter = m->begin(0);
        while ((vert = m->iterate(vertIter)))
        {
          bool allTagged = true;
          m->getIntTag(vert, tag, vals);
          if (vals[0] == 0)
            allTagged = false;

          if (allTagged)
          {
            std::cout << "Flagged vertices " << 1 << " " << allTagged << std::endl;
            dval[0] = 1.0;
            apf::setComponents(field, vert, 0, dval);
          }
        }
        auto iter = cgnsBCMap["Vertex"];
        iter.push_back(std::make_pair(fieldName, field));
      }
      else if (cgnsLocation == "EdgeCenter")
      {
        apf::Field *field = nullptr;
        const std::string fieldName = "EdgeBC_" + bcName;
        field = apf::createField(m, fieldName.c_str(), apf::SCALAR, apf::getConstant(1));

        double dval[1];
        dval[0] = 0.0;
        apf::MeshIterator *edgeIter = m->begin(1);
        apf::MeshEntity *edge = nullptr;
        while ((edge = m->iterate(edgeIter)))
        {
          apf::setComponents(field, edge, 0, dval);
        }
        m->end(edgeIter);

        apf::Downward verts;
        int vals[1];
        edgeIter = m->begin(1);
        while ((edge = m->iterate(edgeIter)))
        {
          const auto numVerts = m->getDownward(edge, 0, verts);
          bool allTagged = true;
          for (int i = 0; i < numVerts; i++)
          {
            m->getIntTag(verts[i], tag, vals);
            if (vals[0] == 0)
              allTagged = false;
          }
          if (allTagged)
          {
            //std::cout << "Flagged edges " << numVerts << " " << allTagged << std::endl;
            dval[0] = 1.0;
            apf::setComponents(field, edge, 0, dval);
          }
        }
        auto iter = cgnsBCMap["EdgeCenter"];
        iter.push_back(std::make_pair(fieldName, field));
      }
      else if (cgnsLocation == "FaceCenter")
      {
        PCU_ALWAYS_ASSERT_VERBOSE(true == false, "Can't have a FaceCenter BC in a 2D mesh");
      }
      else if (cgnsLocation == "CellCenter")
      {
        apf::Field *field = nullptr;
        const std::string fieldName = "CellBC_" + bcName;
        field = apf::createField(m, fieldName.c_str(), apf::SCALAR, apf::getConstant(2));

        double dval[1];
        dval[0] = 0.0;
        apf::MeshIterator *cellIter = m->begin(2);
        apf::MeshEntity *cell = nullptr;
        while ((cell = m->iterate(cellIter)))
        {
          apf::setComponents(field, cell, 0, dval);
        }
        m->end(cellIter);

        apf::Downward verts;
        cellIter = m->begin(2);
        int vals[1];
        while ((cell = m->iterate(cellIter)))
        {
          const auto numVerts = m->getDownward(cell, 0, verts);
          bool allTagged = true;
          for (int i = 0; i < numVerts; i++)
          {
            m->getIntTag(verts[i], tag, vals);
            if (vals[0] == 0)
              allTagged = false;
          }
          if (allTagged)
          {
            //std::cout << "Flagged cells " << numVerts << " " << allTagged << std::endl;
            dval[0] = 1.0;
            apf::setComponents(field, cell, 0, dval);
          }
        }
        auto iter = cgnsBCMap["CellCenter"];
        iter.push_back(std::make_pair(fieldName, field));
      }
    }
    else if (m->getDimension() == 1) // working with a 1D mesh
    {
      if (cgnsLocation == "Vertex")
      {
        apf::Field *field = nullptr;
        const std::string fieldName = "VertBC_" + bcName;
        field = apf::createFieldOn(m, fieldName.c_str(), apf::SCALAR);

        double dval[1];
        dval[0] = 0.0;
        apf::MeshIterator *vertIter = m->begin(0);
        apf::MeshEntity *vert = nullptr;
        while ((vert = m->iterate(vertIter)))
        {
          apf::setComponents(field, vert, 0, dval);
        }
        m->end(vertIter);

        int vals[1];
        vertIter = m->begin(0);
        while ((vert = m->iterate(vertIter)))
        {
          bool allTagged = true;
          m->getIntTag(vert, tag, vals);
          if (vals[0] == 0)
            allTagged = false;

          if (allTagged)
          {
            std::cout << "Flagged vertices " << 1 << " " << allTagged << std::endl;
            dval[0] = 1.0;
            apf::setComponents(field, vert, 0, dval);
          }
        }
        auto iter = cgnsBCMap["Vertex"];
        iter.push_back(std::make_pair(fieldName, field));
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
        apf::Field *field = nullptr;
        const std::string fieldName = "CellBC_" + bcName;
        field = apf::createField(m, fieldName.c_str(), apf::SCALAR, apf::getConstant(1));

        double dval[1];
        dval[0] = 0.0;
        apf::MeshIterator *cellIter = m->begin(1);
        apf::MeshEntity *cell = nullptr;
        while ((cell = m->iterate(cellIter)))
        {
          apf::setComponents(field, cell, 0, dval);
        }
        m->end(cellIter);

        apf::Downward verts;
        cellIter = m->begin(1);
        int vals[1];
        while ((cell = m->iterate(cellIter)))
        {
          const auto numVerts = m->getDownward(cell, 0, verts);
          bool allTagged = true;
          for (int i = 0; i < numVerts; i++)
          {
            m->getIntTag(verts[i], tag, vals);
            if (vals[0] == 0)
              allTagged = false;
          }
          if (allTagged)
          {
            //std::cout << "Flagged cells " << numVerts << " " << allTagged << std::endl;
            dval[0] = 1.0;
            apf::setComponents(field, cell, 0, dval);
          }
        }
        auto iter = cgnsBCMap["CellCenter"];
        iter.push_back(std::make_pair(fieldName, field));
      }
    }

    if (!debugOutput)
      Clean(m);
  }
}; // namespace

void ReadBCInfo(const int cgid, const int base, const int zone, const int nBocos, const int physDim, const int cellDim, const int nsections, std::vector<BCInfo> &bcInfos, const apf::GlobalToVert &globalToVert)
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
      Kill(cgid, "Failed cg_boco_info");

    if (bcMeta.ptsetType == CGNS_ENUMV(PointList) || (bcMeta.ptsetType == CGNS_ENUMV(PointRange)))
    {
      bcMeta.bocoName = std::string(bcMeta.bocoName.c_str());
      //boost::algorithm::trim(bcMeta.bocoName); // can't be bothered including boost

      if (cg_boco_gridlocation_read(cgid, base, zone, boco, &bcMeta.location))
        Kill(cgid, "Failed cg_boco_gridlocation_read");

      bcMeta.locationName = cg_GridLocationName(bcMeta.location);

      if (bcMeta.ptsetType == CGNS_ENUMV(PointRange))
        pointRange = true;
    }
    else
    {
      Kill(cgid,
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
          Kill(cgid, "Failed cg_boco_read");
      }
      else if (bcMeta.locationName == "EdgeCenter") // && cellDim == 2)
      {
        if (cg_boco_read(cgid, base, zone, boco, bcMeta.bcElementIds.data(), NULL))
          Kill(cgid, "Failed cg_boco_read");
      }
      else if (bcMeta.locationName == "FaceCenter") // && cellDim == 3)
      {
        if (cg_boco_read(cgid, base, zone, boco, bcMeta.bcElementIds.data(), NULL))
          Kill(cgid, "Failed cg_boco_read");
      }
      else if (bcMeta.locationName == "CellCenter") // && cellDim == 3)
      {
        if (cg_boco_read(cgid, base, zone, boco, bcMeta.bcElementIds.data(), NULL))
          Kill(cgid, "Failed cg_boco_read");
      }
      else
        Kill(cgid, "Failed Location test for BC Type", bcMeta.locationName,
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
        cgsize_t numElements = -1;
        int verticesPerElement = -1;

        cg_section_read(cgid, base, zone, section, &sectionName[0],
                        &elementType, &el_start, &el_end,
                        &num_bndry, &parent_flag);

        numElements = el_end - el_start + 1;

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

apf::Mesh2 *DoIt(gmi_model *g, const std::string &fname, apf::CGNSBCMap &cgnsBCMap)
{
  int cgid = -1;
  auto comm = PCU_Get_Comm();
  cgp_mpi_comm(comm);
  cgp_pio_mode(CGNS_ENUMV(CGP_INDEPENDENT));
  cgp_open(fname.c_str(), CGNS_ENUMV(CG_MODE_READ), &cgid);

  int nbases = -1;
  cg_nbases(cgid, &nbases);
  if (nbases > 1)
  {
    std::cout << __LINE__ << " CGNS is dead " << std::endl;
    Kill(cgid);
  }

  std::string basename;
  basename.resize(CGIO_MAX_NAME_LENGTH + 1, ' ');

  int cellDim = -1;
  int physDim = -1;
  const int base = 1;
  cg_base_read(cgid, base, &basename[0], &cellDim, &physDim);
  const int readDim = cellDim;

  // Salome cgns is a bit on the odd side: cellDim, physDim, ncoords are not always consistent
  apf::Mesh2 *mesh = apf::makeEmptyMdsMesh(g, cellDim, false);
  apf::GlobalToVert globalToVert;

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
      Kill(cgid);
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
      Kill(cgid);
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

      const auto readElementsAndVerts = [&](const apf::Mesh::Type &type) {
        const auto &ret = ReadElements(cgid, base, zone, section, el_start, el_end, numElements, verticesPerElement);
        if (std::get<1>(ret) > 0)
        {
          const std::vector<cgsize_t> vertexIDs = std::get<0>(ret);
          //apf::construct(mesh, vertexIDs.data(), std::get<1>(ret), type, globalToVert);
          apf::assemble(mesh, vertexIDs.data(), std::get<1>(ret), type, globalToVert); // corresponding finalize below

          const auto nverts = sizes[0];
          const auto ordinates = ReadCGNSCoords(cgid, base, zone, ncoords, nverts, vertexIDs, globalToVert);

          for (const auto &p : globalToVert)
          {
            const auto pp = ordinates.at(p.first);
            apf::Vector3 point(pp[0], pp[1], pp[2]);
            auto iter = globalToVert.find(p.first);
            if (iter != globalToVert.end())
            {
              mesh->setPoint(iter->second, 0, point);
            }
            else
            {
              Kill(cgid, "GlobalToVert lookup problem");
            }
          }
        }
      };

      if (elementType == CGNS_ENUMV(BAR_2))
      {
        if (readDim == 1)
          readElementsAndVerts(apf::Mesh::EDGE);
      }
      else if (elementType == CGNS_ENUMV(QUAD_4))
      {
        if (readDim == 2)
          readElementsAndVerts(apf::Mesh::QUAD);
      }
      else if (elementType == CGNS_ENUMV(TRI_3))
      {
        if (readDim == 2)
          readElementsAndVerts(apf::Mesh::TRIANGLE);
      }
      else if (elementType == CGNS_ENUMV(TETRA_4))
      {
        if (readDim == 3)
          readElementsAndVerts(apf::Mesh::TET);
      }
      else if (elementType == CGNS_ENUMV(PYRA_5))
      {
        if (readDim == 3)
          readElementsAndVerts(apf::Mesh::PYRAMID);
      }
      else if (elementType == CGNS_ENUMV(HEXA_8))
      {
        if (readDim == 3)
          readElementsAndVerts(apf::Mesh::HEX);
      }
      else
      {
        std::cout << __LINE__ << " CGNS is dead "
                  << " " << SupportedCGNSElementTypeToString(elementType) << std::endl;
        Kill(cgid);
      }
    }

    if (nBocos > 0)
    {
      std::cout << "Attempting to read BCS info "
                << " " << nBocos << std::endl;
      ReadBCInfo(cgid, base, zone, nBocos, physDim, cellDim, nsections, bcInfos, globalToVert);
    }
  }

  apf::finalise(mesh, globalToVert);
  apf::alignMdsRemotes(mesh);
  apf::deriveMdsModel(mesh);
  mesh->acceptChanges();
  apf::verify(mesh, true);
  {
    apf::GlobalNumbering *gn = nullptr;
    gn = apf::makeGlobal(apf::numberOwnedNodes(mesh, "vertex Indices"));
    apf::synchronize(gn);
  }
  // no synchronize call
  // https://github.com/SNLComputation/Albany/blob/master/src/disc/pumi/Albany_APFDiscretization.cpp @ various place throughout file
  // https://github.com/SCOREC/core/issues/249
  {
    apf::GlobalNumbering *gn = nullptr;
    gn = apf::makeGlobal(apf::numberElements(mesh, "element Indices"));
  }

  for (auto &bc : bcInfos)
  {
    bc.TagVertices(cgid, mesh, globalToVert);
  }

  for (auto &bc : bcInfos)
  {
    bc.TagBCEntities(cgid, mesh, cgnsBCMap);
  }

  if (PCU_Comm_Initialized())
    cgp_close(cgid);
  else
    cg_close(cgid);

  return mesh;
}
} // namespace

namespace apf
{

// caller needs to bring up and pull down mpi/pcu: mpi/pcu is required and assumed.
Mesh2 *loadMdsFromCGNS(gmi_model *g, const char *fname, apf::CGNSBCMap &cgnsBCMap)
{
#ifdef HAVE_CGNS
  Mesh2 *m = DoIt(g, fname, cgnsBCMap);
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
