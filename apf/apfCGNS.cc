/*
 * Author: Dr Andrew Parker (2019) - FGE Ltd
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfMesh.h"
#include "apfNumbering.h"
#include "apfNumberingClass.h"
#include "apfShape.h"
#include "apfFieldData.h"
#include <pcu_util.h>
#include <lionPrint.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <stdint.h>
#include <array>
#include <vector>
#include <numeric>

// === includes for safe_mkdir ===
#include <reel.h>
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h>  /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h>     /* for checking the error from mkdir */
// ===============================

#ifdef HAVE_CGNS
//
#include <cgns_io.h>
#include <pcgnslib.h>
//
#endif

// Note: currently, even in 2D or 1D a full 3D vector and matrix are written out to the file

namespace
{
//https://proteustoolkit.org/capi/html/_parallel_mesh_converter_8cpp_source.html
using Count = std::pair<int, int>;
static Count count(apf::Mesh *m, int dim)
{
  const int local = apf::countOwned(m, dim);
  int total = local;
  m->getPCU()->Add<int>(&total, 1); // size of total array
  return std::make_pair(total, local);
}

struct CGNS
{
  int index = -1;
  int base = -1;
  int zone = -1;
  const int phys_dim = 3;
  std::string fname;
};

/*
* Disable this until I find out why all tags get an _tet_ or _pyr_ fieldname, makes the output file super messy, and fields seem duplicated in tags
void WriteTags(const CGNS &cgns, const std::vector<std::vector<apf::MeshEntity *>> &orderedEnts, const std::vector<std::pair<cgsize_t, cgsize_t>> &ranges, const std::vector<apf::MeshEntity *> &orderedVertices, const int &vStart, const int &vEnd, apf::Mesh *m)
{
  const auto loopVertexTags = [&m](const auto &orderedEnts, const int &solIndex, const auto &inner, const auto &post, const int &start, const int &end) {
    apf::DynamicArray<apf::MeshTag *> tags;
    m->getTags(tags);
    for (std::size_t i = 0; i < tags.getSize(); ++i)
    {
      apf::MeshTag *t = tags[i];
      const int tagType = m->getTagType(t);
      const int tagSize = m->getTagSize(t);
      std::string tagName(m->getTagName(t));

      if (tagName.size() > 32)
        tagName.resize(32); // daft api

      // boring... replace with variant
      std::vector<int> idata;
      std::vector<double> ddata;
      std::vector<long> ldata;

      for (int ti = 0; ti < tagSize; ti++)
      {
        cgsize_t rmin[3];
        cgsize_t rmax[3];

        rmin[0] = start;
        rmax[0] = end;
        //std::cout << start << " " << end << std::endl;
        int fieldIndex = -1;
        auto tagNameNew = tagName;
        if (tagSize > 1)
          tagNameNew += "_[" + std::to_string(ti) + "]";

        for (const auto &e : orderedEnts)
        {
          if (m->hasTag(e, t) && m->isOwned(e))
          {
            inner(e, t, tagSize, ti, idata, ddata, ldata, tagType);
          }
        }

        // ensure collectives are called by all, even if local has no data
        int isize = idata.size();
        PCU_Add_Ints(&isize, 1); // size of total array

        int dsize = ddata.size();
        PCU_Add_Ints(&dsize, 1); // size of total array

        int lsize = ldata.size();
        PCU_Add_Ints(&lsize, 1); // size of total array

        // oddness of the api
        rmin[1] = rmin[0];
        rmin[2] = rmin[0];
        rmax[1] = rmax[0];
        rmax[2] = rmax[0];

        if (isize > 0 || dsize > 0 || lsize > 0)
          post(solIndex, tagNameNew, idata, ddata, ldata, rmin, rmax, isize, dsize, lsize, fieldIndex);
      }
    }
  };

  const auto loopCellTags = [&m](const auto &orderedEnts, const int &solIndex, const auto &inner, const auto &post, const auto &ranges) {
    apf::DynamicArray<apf::MeshTag *> tags;
    m->getTags(tags);
    for (std::size_t i = 0; i < tags.getSize(); ++i)
    {
      apf::MeshTag *t = tags[i];
      const int tagType = m->getTagType(t);
      const int tagSize = m->getTagSize(t);
      std::string tagName(m->getTagName(t));

      if (tagName.size() > 32)
        tagName.resize(32); // daft api

      for (int ti = 0; ti < tagSize; ti++)
      {
        int fieldIndex = -1;
        auto tagNameNew = tagName;
        if (tagSize > 1)
          tagNameNew += "_[" + std::to_string(ti) + "]";

        for (std::size_t e = 0; e < orderedEnts.size(); e++)
        {
          // boring... replace with variant
          std::vector<int> idata;
          std::vector<double> ddata;
          std::vector<long> ldata;

          cgsize_t rmin[3];
          cgsize_t rmax[3];

          rmin[0] = ranges[e].first;
          rmax[0] = ranges[e].second;

          for (const auto &elm : orderedEnts[e])
          {
            if (m->hasTag(elm, t) && m->isOwned(elm))
            {
              inner(elm, t, tagSize, ti, idata, ddata, ldata, tagType);
            }
          }

          // ensure collectives are called by all, even if local has no data
          int isize = idata.size();
          PCU_Add_Ints(&isize, 1); // size of total array

          int dsize = ddata.size();
          PCU_Add_Ints(&dsize, 1); // size of total array

          int lsize = ldata.size();
          PCU_Add_Ints(&lsize, 1); // size of total array

          // oddness of the api
          rmin[1] = rmin[0];
          rmin[2] = rmin[0];
          rmax[1] = rmax[0];
          rmax[2] = rmax[0];

          if (isize > 0 || dsize > 0 || lsize > 0)
            post(solIndex, tagNameNew, idata, ddata, ldata, rmin, rmax, isize, dsize, lsize, fieldIndex);
        }
      }
    }
  };

  const auto postLambda = [&cgns](const int &solIndex, const std::string &name, std::vector<int> &idata, std::vector<double> &ddata, std::vector<long> &ldata, const cgsize_t *rmin, const cgsize_t *rmax, const int isize, const int dsize, const int lsize, int &fieldIndex) {
    if (dsize > 0)
    {
      if (fieldIndex == -1)
      {
        if (cgp_field_write(cgns.index, cgns.base, cgns.zone, solIndex, CGNS_ENUMV(RealDouble), name.c_str(), &fieldIndex))
          cgp_error_exit();
      }
      if (cgp_field_write_data(cgns.index, cgns.base, cgns.zone, solIndex, fieldIndex, &rmin[0], &rmax[0],
                               ddata.data()))
        cgp_error_exit();
    }
    else if (isize > 0)
    {
      if (fieldIndex == -1)
      {
        if (cgp_field_write(cgns.index, cgns.base, cgns.zone, solIndex, CGNS_ENUMV(Integer), name.c_str(), &fieldIndex))
          cgp_error_exit();
      }

      if (cgp_field_write_data(cgns.index, cgns.base, cgns.zone, solIndex, fieldIndex, &rmin[0], &rmax[0],
                               idata.data()))
        cgp_error_exit();
    }
    else if (lsize > 0)
    {
      if (fieldIndex == -1)
      {
        if (cgp_field_write(cgns.index, cgns.base, cgns.zone, solIndex, CGNS_ENUMV(LongInteger), name.c_str(), &fieldIndex))
          cgp_error_exit();
      }
      if (cgp_field_write_data(cgns.index, cgns.base, cgns.zone, solIndex, fieldIndex, &rmin[0], &rmax[0],
                               ldata.data()))
        cgp_error_exit();
    }
  };

  const auto innerLambda = [&m](apf::MeshEntity *elem, apf::MeshTag *tag, const int &tagSize, const int &ti, std::vector<int> &idata, std::vector<double> &ddata, std::vector<long> &ldata, const int &tagType) {
    if (tagType == apf::Mesh::TagType::DOUBLE)
    {
      std::vector<double> vals(tagSize, -12345);
      m->getDoubleTag(elem, tag, vals.data());
      ddata.push_back(vals[ti]);
    }
    else if (tagType == apf::Mesh::TagType::INT)
    {
      std::vector<int> vals(tagSize, -12345);
      m->getIntTag(elem, tag, vals.data());
      idata.push_back(vals[ti]);
    }
    else if (tagType == apf::Mesh::TagType::LONG)
    {
      std::vector<long> vals(tagSize, -12345);
      m->getLongTag(elem, tag, vals.data());
      ldata.push_back(vals[ti]);
    }
    else
    {
      std::cout << "Strange" << std::endl;
      exit(-1);
    }
  };

  int solIndex = -1;

  {
    if (cg_sol_write(cgns.index, cgns.base, cgns.zone, "Vertex Tag Data", CGNS_ENUMV(Vertex), &solIndex))
      cg_error_exit();

    loopVertexTags(orderedVertices, solIndex, innerLambda, postLambda, vStart, vEnd);
  }

  {
    if (cg_sol_write(cgns.index, cgns.base, cgns.zone, "Cell Tag Data", CGNS_ENUMV(CellCenter), &solIndex))
      cg_error_exit();

    loopCellTags(orderedEnts, solIndex, innerLambda, postLambda, ranges);
  }
}
*/
void WriteFields(const CGNS &cgns, const std::vector<std::vector<apf::MeshEntity *>> &orderedEnts, const std::vector<std::pair<cgsize_t, cgsize_t>> &ranges, const std::vector<apf::MeshEntity *> &orderedVertices, const int &vStart, const int &vEnd, apf::Mesh *m)
{
  const auto writeField = [&m, &cgns](apf::Field *f, const auto &orderedEnts, const int &solIndex, const auto &inner, const auto &post, const int &numComponents, const int &component, const std::string &fieldName, const int &start, const int &end, int &fieldIndex) {
    std::vector<double> data;

    cgsize_t rmin[3];
    cgsize_t rmax[3];

    rmin[0] = start;
    rmax[0] = end;

    apf::FieldDataOf<double> *fieldData = f->getData();
    for (const auto &e : orderedEnts)
    {
      if (fieldData->hasEntity(e) && m->isOwned(e))
      {
        inner(e, fieldData, data, numComponents, component);
      }
    }

    int size = data.size();
    m->getPCU()->Add<int>(&size, 1); // size of total array

    // oddness of the api
    rmin[1] = rmin[0];
    rmin[2] = rmin[0];
    rmax[1] = rmax[0];
    rmax[2] = rmax[0];
    if (size > 0)
    {
      if (fieldIndex == -1)
      {
        //std::cout << "FieldName " << fieldName << std::endl;
        if (cgp_field_write(cgns.index, cgns.base, cgns.zone, solIndex, CGNS_ENUMV(RealDouble), fieldName.c_str(), &fieldIndex))
          cgp_error_exit();
      }

      post(solIndex, data, rmin, rmax, size, fieldIndex);
    }
  };

  const auto loopCellFields = [&m, &writeField](const auto &orderedEnts, const int &solIndex, const auto &inner, const auto &post, const auto &ranges) {
    for (int i = 0; i < m->countFields(); ++i)
    {
      apf::Field *f = m->getField(i);
      const int numComponents = f->countComponents();
      std::string fieldName(f->getName());
      if (fieldName.size() > 32)
        fieldName.resize(32); // daft api

      for (int component = 0; component < numComponents; component++)
      {
        int fieldIndex = -1;
        auto fieldNameNew = fieldName;
        if (numComponents > 1)
          fieldNameNew += "_[" + std::to_string(component) + "]";

        for (std::size_t e = 0; e < orderedEnts.size(); e++)
        {
          //std::cout << "CELL " << fieldNameNew << " " << fieldName << " " << component << " " << numComponents << " " << e << " " << ranges[e].first << " " << ranges[e].second << std::endl;
          writeField(f, orderedEnts[e], solIndex, inner, post, numComponents, component, fieldNameNew, ranges[e].first, ranges[e].second, fieldIndex);
        }
      }
    }
  };

  const auto loopVertexFields = [&m, &writeField](const auto &orderedEnts, const int &solIndex, const auto &inner, const auto &post, const int &vStart, const int &vEnd) {
    for (int i = 0; i < m->countFields(); ++i)
    {
      apf::Field *f = m->getField(i);
      const int numComponents = f->countComponents();
      std::string fieldName(f->getName());
      if (fieldName.size() > 32)
        fieldName.resize(32); // daft api

      for (int component = 0; component < numComponents; component++)
      {
        int fieldIndex = -1;
        auto fieldNameNew = fieldName;
        if (numComponents > 1)
          fieldNameNew += "_[" + std::to_string(component) + "]";

        //std::cout << "VERTEX " << fieldNameNew << " " << fieldName << " " << component << " " << numComponents << std::endl;
        writeField(f, orderedEnts, solIndex, inner, post, numComponents, component, fieldNameNew, vStart, vEnd, fieldIndex);
      }
    }
  };

  const auto postLambda = [&cgns](const int &solIndex, std::vector<double> &ddata, const cgsize_t *rmin, const cgsize_t *rmax, const int &globalSize, const int &fieldIndex) {
    if (globalSize > 0)
    {
      if (cgp_field_write_data(cgns.index, cgns.base, cgns.zone, solIndex, fieldIndex, &rmin[0], &rmax[0],
                               ddata.data()))
        cgp_error_exit();
    }
  };

  const auto innerLambda = [](apf::MeshEntity *elem, apf::FieldDataOf<double> *fieldData, std::vector<double> &ddata, const int &numComponents, const int &component) {
    std::vector<double> vals(numComponents, -12345);
    fieldData->get(elem, vals.data());
    //std::cout << numComponents << " " << component << " " << vals[0] << std::endl;
    ddata.push_back(vals[component]);
  };

  int solIndex = -1;

  {
    if (cg_sol_write(cgns.index, cgns.base, cgns.zone, "Vertex Field Data", CGNS_ENUMV(Vertex), &solIndex))
      cg_error_exit();

    loopVertexFields(orderedVertices, solIndex, innerLambda, postLambda, vStart, vEnd);
  }

  {
    if (cg_sol_write(cgns.index, cgns.base, cgns.zone, "Cell Field Data", CGNS_ENUMV(CellCenter), &solIndex))
      cg_error_exit();

    loopCellFields(orderedEnts, solIndex, innerLambda, postLambda, ranges);
  }
}

auto WriteVertices(const CGNS &cgns, apf::Mesh *m, apf::GlobalNumbering *gvn)
{
  int Cx = -1;
  int Cy = -1;
  int Cz = -1;

  if (cgns.phys_dim > 0)
  {
    if (cgp_coord_write(cgns.index, cgns.base, cgns.zone, CGNS_ENUMV(RealDouble), "CoordinateX", &Cx))
      cgp_error_exit();
  }
  if (cgns.phys_dim > 1)
  {
    if (cgp_coord_write(cgns.index, cgns.base, cgns.zone, CGNS_ENUMV(RealDouble), "CoordinateY", &Cy))
      cgp_error_exit();
  }
  if (cgns.phys_dim > 2)
  {
    if (cgp_coord_write(cgns.index, cgns.base, cgns.zone, CGNS_ENUMV(RealDouble), "CoordinateZ", &Cz))
      cgp_error_exit();
  }

  std::vector<apf::MeshEntity *> orderedVertices;
  cgsize_t vertexMin[3];
  cgsize_t vertexMax[3];
  cgsize_t contigRange = -1;
  {
    std::array<std::vector<double>, 3> coords;

    vertexMin[0] = std::numeric_limits<cgsize_t>::max();
    vertexMax[0] = 0;

    apf::Vector3 point;
    apf::MeshIterator *vertIter = m->begin(0);
    apf::MeshEntity *vert = nullptr;
    while ((vert = m->iterate(vertIter)))
    {
      if (m->isOwned(vert))
      {
        const cgsize_t n = static_cast<cgsize_t>(apf::getNumber(gvn, vert, 0) + 1); // one based
        if (contigRange == -1)
        {
          contigRange = n;
        }
        else
        {
          const auto predict = contigRange + 1;

          if (n != predict)
          {
            //std::cout << std::string("Range must be contigious for vertices ") + std::to_string(n) + " " + std::to_string(contigRange) << std::endl;
            PCU_ALWAYS_ASSERT_VERBOSE(true == false, (std::string("Range must be contigious for vertices ") + std::to_string(n) + " " + std::to_string(contigRange)).c_str());
          }
          else
            contigRange = predict; // == n
        }

        vertexMin[0] = std::min(vertexMin[0], n);
        vertexMax[0] = std::max(vertexMax[0], n);

        m->getPoint(vert, 0, point);
        coords[0].push_back(point[0]);
        coords[1].push_back(point[1]);
        coords[2].push_back(point[2]);

        orderedVertices.push_back(vert);
      }
    }
    m->end(vertIter);

    // oddness of the api
    vertexMin[1] = vertexMin[0];
    vertexMin[2] = vertexMin[0];
    vertexMax[1] = vertexMax[0];
    vertexMax[2] = vertexMax[0];

    if (cgns.phys_dim > 0)
    {
      if (cgp_coord_write_data(cgns.index, cgns.base, cgns.zone, Cx, &vertexMin[0], &vertexMax[0], coords[0].data()))
        cgp_error_exit();
    }
    if (cgns.phys_dim > 1)
    {
      if (cgp_coord_write_data(cgns.index, cgns.base, cgns.zone, Cy, &vertexMin[0], &vertexMax[0], coords[1].data()))
        cgp_error_exit();
    }
    if (cgns.phys_dim > 2)
    {
      if (cgp_coord_write_data(cgns.index, cgns.base, cgns.zone, Cz, &vertexMin[0], &vertexMax[0], coords[2].data()))
        cgp_error_exit();
    }
  }
  return std::make_tuple(orderedVertices, vertexMin[0], vertexMax[0]);
}

using CellElementReturn = std::tuple<std::vector<std::vector<apf::MeshEntity *>>, std::vector<std::pair<cgsize_t, cgsize_t>>>;
CellElementReturn WriteElements(const CGNS &cgns, apf::Mesh *m, apf::GlobalNumbering *gvn, const int cell_dim, const Count &cellCount, const std::vector<apf::Mesh::Type> &apfElementOrder, const std::vector<CGNS_ENUMT(ElementType_t)> &cgnsElementOrder)
{
  std::vector<int> globalNumbersByElementType(apfElementOrder.size(), 0);
  std::vector<int> numbersByElementType(apfElementOrder.size(), 0);
  for (std::size_t o = 0; o < apfElementOrder.size(); o++)
  {
    apf::MeshIterator *cellIter = m->begin(cell_dim);
    apf::MeshEntity *cell = nullptr;
    int counter = 0;
    while ((cell = m->iterate(cellIter)))
    {
      if (m->getType(cell) == apfElementOrder[o] && m->isOwned(cell)) // must be same test as below
      {
        counter++;
      }
    }
    m->end(cellIter);
    numbersByElementType[o] = counter;
    int total = counter;
    m->getPCU()->Add<int>(&total, 1); // size of total array
    globalNumbersByElementType[o] = total;
  }
  cgsize_t allTotal = std::accumulate(globalNumbersByElementType.begin(), globalNumbersByElementType.end(), 0);
  PCU_ALWAYS_ASSERT_VERBOSE(allTotal == cellCount.first, ("Must be equal " + std::to_string(allTotal) + " " + std::to_string(cellCount.first)).c_str());

  int globalStart = 1; // one-based
  std::vector<std::vector<apf::MeshEntity *>> orderedElements(apfElementOrder.size());
  std::vector<std::pair<cgsize_t, cgsize_t>> ranges(apfElementOrder.size());
  for (std::size_t o = 0; o < apfElementOrder.size(); o++)
  {
    std::vector<cgsize_t> elementVertices;
    apf::MeshIterator *cellIter = m->begin(cell_dim);
    apf::MeshEntity *cell = nullptr;
    apf::Downward verts;

    while ((cell = m->iterate(cellIter)))
    {
      if (m->getType(cell) == apfElementOrder[o] && m->isOwned(cell)) // must be same test as above
      {
        const auto numVerts = m->getDownward(cell, 0, verts);
        for (int i = 0; i < numVerts; i++)
        {
          const auto n = apf::getNumber(gvn, verts[i], 0);
          elementVertices.push_back(n + 1); // one-based
        }
        orderedElements[o].push_back(cell);
      }
    }
    m->end(cellIter);

    if (globalNumbersByElementType[o] > 0)
    {
      const int globalEnd = globalStart + globalNumbersByElementType[o] - 1; // one-based stuff
      //
      int sectionNumber = -1;
      const std::string name = std::string(cg_ElementTypeName(cgnsElementOrder[o])) + " " + std::to_string(globalStart) + "->" + std::to_string(globalEnd);
      //std::cout << "Section Name " << name << std::endl;
      if (cgp_section_write(cgns.index, cgns.base, cgns.zone, name.c_str(), cgnsElementOrder[o], globalStart, globalEnd, 0, &sectionNumber)) // global start, end within the file for that element type
        cgp_error_exit();

      std::vector<int> allNumbersForThisType(m->getPCU()->Peers(), 0);
      MPI_Allgather(&numbersByElementType[o], 1, MPI_INT, allNumbersForThisType.data(), 1,
                    MPI_INT, m->getPCU()->GetMPIComm());

      cgsize_t num = 0;
      for (int i = 0; i < m->getPCU()->Self(); i++)
        num += allNumbersForThisType[i];

      cgsize_t elStart = globalStart + num;
      cgsize_t elEnd = elStart + numbersByElementType[o] - 1;                                      // one-based stuff
      if (cgp_elements_write_data(cgns.index, cgns.base, cgns.zone, sectionNumber, elStart, elEnd, // per processor within the range[start, end]
                                  elementVertices.data()))
        cgp_error_exit();

      //std::cout << std::flush << "RANK: " << PCU_Comm_Self() << " ==> " << globalStart << " " << globalEnd << " elStart " << elStart << " elEnd " << elEnd << " numbersByElementType[o] " << numbersByElementType[o] << " of " << std::string(cg_ElementTypeName(cgnsElementOrder[o])) << std::flush << std::endl;

      globalStart += globalNumbersByElementType[o];

      ranges[o] = std::make_pair(elStart, elEnd);
    }
    else
      ranges[o] = std::make_pair(-1, -1);
  }
  return std::make_tuple(orderedElements, ranges);
}

void AddBocosToMainBase(const CGNS &cgns, const CellElementReturn &cellResults, const int &cellCount, apf::Mesh *m, const apf::CGNSBCMap &cgnsBCMap, const std::map<apf::Mesh::Type, CGNS_ENUMT(ElementType_t)> &apf2cgns, apf::GlobalNumbering *gvn)
{
  const auto EdgeLoop = [&m](const auto &lambda, apf::MeshTag *edgeTag) {
    apf::MeshIterator *edgeIter = m->begin(1);
    apf::MeshEntity *edge = nullptr;
    int vals[1];
    vals[0] = -1;

    while ((edge = m->iterate(edgeIter)))
    {
      m->getIntTag(edge, edgeTag, vals);
      if (vals[0] == 1 && m->isOwned(edge))
      {
        lambda(edge);
      }
    }
    m->end(edgeIter);
  };

  const auto FaceLoop = [&m](const auto &lambda, apf::MeshTag *faceTag) {
    apf::MeshIterator *faceIter = m->begin(2);
    apf::MeshEntity *face = nullptr;
    int vals[1];
    vals[0] = -1;

    while ((face = m->iterate(faceIter)))
    {
      m->getIntTag(face, faceTag, vals);
      if (vals[0] == 1 && m->isOwned(face))
      {
        lambda(face);
      }
    }
    m->end(faceIter);
  };

  const auto BCEntityAdder = [&apf2cgns, &m, &cgns, &gvn](const auto &Looper, const auto &bcGroup, int &startingLocation) {
    std::map<apf::Mesh::Type, std::vector<apf::MeshEntity *>> bcEntTypes;
    for (const auto &r : apf2cgns)
      bcEntTypes.insert(std::make_pair(r.first, std::vector<apf::MeshEntity *>()));

    const auto lambda = [&m, &bcEntTypes](apf::MeshEntity *entity) {
      const auto type = m->getType(entity);
      auto iter = bcEntTypes.find(type);
      if (iter != bcEntTypes.end())
      {
        iter->second.push_back(entity);
      }
      else
      {
        PCU_ALWAYS_ASSERT_VERBOSE(true == false, "must not come in here");
      }
    };
    //
    Looper(lambda, bcGroup.bcsMarkerTag);
    //
    const int cacheStart = startingLocation + 1;
    int cacheEnd = -1;
    for (const auto &bc : bcEntTypes)
    {
      int startOfBCBlock = startingLocation + 1;
      const int number = bc.second.size();
      int total = number;
      m->getPCU()->Add<int<(&total, 1); // size of total array
      if (total > 0)
      {
        const auto allEnd = startOfBCBlock + total - 1; //one-based
        int sectionNumber = -1;
        const std::string name = bcGroup.cgnsBCSName + " " + std::to_string(startOfBCBlock) + "->" + std::to_string(allEnd); //cg_ElementTypeName(apf2cgns[bc.first]);
        if (cgp_section_write(cgns.index, cgns.base, cgns.zone, name.c_str(), apf2cgns.at(std::get<0>(bc)), startOfBCBlock, allEnd, 0,
                              &sectionNumber))
          cgp_error_exit();

        std::vector<cgsize_t> elementVertices;
        for (std::size_t e = 0; e < bc.second.size(); e++)
        {
          apf::Downward verts;
          const auto numVerts = m->getDownward(bc.second[e], 0, verts);
          for (int i = 0; i < numVerts; i++)
          {
            const auto n = apf::getNumber(gvn, verts[i], 0);
            elementVertices.push_back(n + 1); // one-based
          }
        }

        std::vector<int> allNumbersForThisType(m->getPCU()->Peers(), 0);
        MPI_Allgather(&number, 1, MPI_INT, allNumbersForThisType.data(), 1,
                      MPI_INT, m->getPCU()->GetMPIComm());

        cgsize_t num = 0;
        for (int i = 0; i < m->getPCU()->Self(); i++)
          num += allNumbersForThisType[i];

        cgsize_t elStart = startOfBCBlock + num;
        cgsize_t elEnd = elStart + number - 1; // one-based stuff

        if (number == 0)
        {
          if (cgp_elements_write_data(cgns.index, cgns.base, cgns.zone, sectionNumber, elStart, elEnd, nullptr))
            cgp_error_exit();
        }
        else
        {
          if (cgp_elements_write_data(cgns.index, cgns.base, cgns.zone, sectionNumber, elStart, elEnd,
                                      elementVertices.data()))
            cgp_error_exit();
        }

        // Not parallel correct
        startingLocation = allEnd;
        cacheEnd = allEnd;
      }
    }
    std::vector<int> cacheStarts(m->getPCU()->Peers(), 0);
    MPI_Allgather(&cacheStart, 1, MPI_INT, cacheStarts.data(), 1,
                  MPI_INT, m->getPCU()->GetMPIComm());
    std::vector<int> cacheEnds(m->getPCU()->Peers(), 0);
    MPI_Allgather(&cacheEnd, 1, MPI_INT, cacheEnds.data(), 1,
                  MPI_INT, m->getPCU()->GetMPIComm());
    return std::make_pair(cacheStarts, cacheEnds);
  };

  const auto globalElementList = [](const std::vector<cgsize_t> &bcList, std::vector<cgsize_t> &allElements) {
    std::vector<int> sizes(m->getPCU()->Peers(), 0); // important initialiser
    const int l = bcList.size();
    MPI_Allgather(&l, 1, MPI_INT, sizes.data(), 1,
                  MPI_INT, m->getPCU()->GetMPIComm());

    int totalLength = 0;
    for (const auto &i : sizes)
      totalLength += i;

    std::vector<int> displacement(m->getPCU()->Peers(), -1);
    displacement[0] = 0;
    for (std::size_t i = 1; i < displacement.size(); i++)
      displacement[i] = displacement[i - 1] + sizes[i - 1];

    allElements.resize(totalLength);
    MPI_Allgatherv(bcList.data(), bcList.size(), MPI_INT, allElements.data(),
                   sizes.data(), displacement.data(), MPI_INT,
                   m->getPCU()->GetMPIComm());
  };

  const auto doVertexBC = [&](const auto &iter) {
    for (const auto &p : iter->second)
    {
      std::vector<cgsize_t> bcList;
      apf::MeshIterator *vertIter = m->begin(0);
      apf::MeshEntity *vert = nullptr;
      int vals[1];
      vals[0] = -1;

      while ((vert = m->iterate(vertIter)))
      {
        m->getIntTag(vert, p.bcsMarkerTag, vals);
        if (vals[0] == 1 && m->isOwned(vert))
        {
          const auto n = apf::getNumber(gvn, vert, 0);
          bcList.push_back(n + 1); // one-based
        }
      }
      m->end(vertIter);

      // Annoying that the API requires this!
      std::vector<cgsize_t> allElements;
      globalElementList(bcList, allElements);

      int bcIndex = -1;
      if (cg_boco_write(cgns.index, cgns.base, cgns.zone, p.cgnsBCSName.c_str(), CGNS_ENUMV(BCGeneral), CGNS_ENUMV(PointList), allElements.size(),
                        allElements.data(), &bcIndex))
        cg_error_exit();

      CGNS_ENUMT(GridLocation_t)
      location = CGNS_ENUMV(Vertex);
      if (cg_boco_gridlocation_write(cgns.index, cgns.base, cgns.zone, bcIndex, location))
        cg_error_exit();
    }
  };

  const auto doEdgeBC = [&](const auto &iter, int &startingLocation) {
    for (const auto &p : iter->second)
    {
      const auto se = BCEntityAdder(EdgeLoop, p, startingLocation);
      for (int i = 0; i < m->getPCU()->Peers(); i++)
      {
        PCU_ALWAYS_ASSERT_VERBOSE(se.first[i] == se.first[0], "Must all be the same ");
        PCU_ALWAYS_ASSERT_VERBOSE(se.second[i] == se.second[0], "Must all be the same ");
      }

      const std::array<cgsize_t, 2> bcRange = {{se.first[0], se.second[0]}};
      int bcIndex = -1;
      if (cg_boco_write(cgns.index, cgns.base, cgns.zone, p.cgnsBCSName.c_str(), CGNS_ENUMV(BCGeneral), CGNS_ENUMV(PointRange), 2,
                        bcRange.data(), &bcIndex))
        cg_error_exit();

      CGNS_ENUMT(GridLocation_t)
      location = CGNS_ENUMV(EdgeCenter);
      if (cg_boco_gridlocation_write(cgns.index, cgns.base, cgns.zone, bcIndex, location))
        cg_error_exit();
    }
  };

  const auto doFaceBC = [&](const auto &iter, int &startingLocation) {
    for (const auto &p : iter->second)
    {
      const auto se = BCEntityAdder(FaceLoop, p, startingLocation);

      for (int i = 0; i < m->getPCU()->Peers(); i++)
      {
        PCU_ALWAYS_ASSERT_VERBOSE(se.first[i] == se.first[0], "Must all be the same ");
        PCU_ALWAYS_ASSERT_VERBOSE(se.second[i] == se.second[0], "Must all be the same ");
      }
      const std::array<cgsize_t, 2> bcRange = {{se.first[0], se.second[0]}};
      int bcIndex = -1;
      if (cg_boco_write(cgns.index, cgns.base, cgns.zone, p.cgnsBCSName.c_str(), CGNS_ENUMV(BCGeneral), CGNS_ENUMV(PointRange), 2,
                        bcRange.data(), &bcIndex))
        cg_error_exit();

      CGNS_ENUMT(GridLocation_t)
      location = CGNS_ENUMV(FaceCenter);
      if (cg_boco_gridlocation_write(cgns.index, cgns.base, cgns.zone, bcIndex, location))
        cg_error_exit();
    }
  };

  const auto doCellBC = [&](const auto &iter, const int &) {
    for (const auto &p : iter->second)
    {
      std::vector<cgsize_t> bcList;
      int vals[1];
      vals[0] = -1;

      const auto &elements = std::get<0>(cellResults);
      const auto &ranges = std::get<1>(cellResults);

      for (std::size_t e = 0; e < elements.size(); e++)
      {
        std::size_t counter = ranges[e].first;
        for (const auto &elm : elements[e])
        {
          m->getIntTag(elm, p.bcsMarkerTag, vals);
          if (vals[0] == 1 && m->isOwned(elm))
          {
            bcList.push_back(counter);
          }
          counter++;
        }
      }

      // Annoying that the API requires this!
      std::vector<cgsize_t> allElements;
      globalElementList(bcList, allElements);

      int bcIndex = -1;
      if (cg_boco_write(cgns.index, cgns.base, cgns.zone, p.cgnsBCSName.c_str(), CGNS_ENUMV(BCGeneral), CGNS_ENUMV(PointList), allElements.size(),
                        allElements.data(), &bcIndex))
        cg_error_exit();

      CGNS_ENUMT(GridLocation_t)
      location = CGNS_ENUMV(CellCenter);
      if (cg_boco_gridlocation_write(cgns.index, cgns.base, cgns.zone, bcIndex, location))
        cg_error_exit();
    }
  };

  const auto cell_dim = m->getDimension();
  if (cell_dim == 3)
  {
    int startingLocation = cellCount;

    auto iter = cgnsBCMap.find("Vertex");
    if (iter != cgnsBCMap.end())
    {
      doVertexBC(iter);
    }
    iter = cgnsBCMap.find("EdgeCenter");
    if (iter != cgnsBCMap.end())
    {
      doEdgeBC(iter, startingLocation);
    }
    iter = cgnsBCMap.find("FaceCenter");
    if (iter != cgnsBCMap.end())
    {
      doFaceBC(iter, startingLocation);
    }
    iter = cgnsBCMap.find("CellCenter");
    if (iter != cgnsBCMap.end())
    {
      doCellBC(iter, 3);
    }
  }
  else if (cell_dim == 2)
  {
    int startingLocation = cellCount;

    auto iter = cgnsBCMap.find("Vertex");
    if (iter != cgnsBCMap.end())
    {
      doVertexBC(iter);
    }
    iter = cgnsBCMap.find("EdgeCenter");
    if (iter != cgnsBCMap.end())
    {
      doEdgeBC(iter, startingLocation);
    }
    iter = cgnsBCMap.find("CellCenter");
    if (iter != cgnsBCMap.end())
    {
      doCellBC(iter, 2);
    }
  }
  else if (cell_dim == 1)
  {
    auto iter = cgnsBCMap.find("Vertex");
    if (iter != cgnsBCMap.end())
    {
      doVertexBC(iter);
    }
    iter = cgnsBCMap.find("CellCenter");
    if (iter != cgnsBCMap.end())
    {
      doCellBC(iter, 1);
    }
  }
}

// copy is deliberate
void Write3DFaces(CGNS cgns, apf::Mesh *m, const Count &faceCount, const Count &vertexCount, apf::GlobalNumbering *gvn)
{
  std::array<cgsize_t, 3> sizes;
  sizes[0] = vertexCount.first; // global
  sizes[1] = faceCount.first;   // global
  sizes[2] = 0;                 // nodes are unsorted, as defined by api

  {
    std::string baseName("FaceMeshBase");
    if (cg_base_write(cgns.index, baseName.c_str(), 2, cgns.phys_dim, &cgns.base))
      cg_error_exit();
  }
  // Write the default units at the top.
  if (cg_goto(cgns.index, cgns.base, "end"))
    cg_error_exit();

  if (cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin),
                     CGNS_ENUMV(Degree)))
    cg_error_exit();

  if (cg_dataclass_write(CGNS_ENUMV(Dimensional)))
    cg_error_exit();

  {
    std::string zoneName("FaceMeshZone");
    if (cg_zone_write(cgns.index, cgns.base, zoneName.c_str(), sizes.data(), CGNS_ENUMV(Unstructured), &cgns.zone))
      cg_error_exit();
  }

  // Stupid api why in the name would you have to repeat ALL vertices again....
  auto &&vertResult = WriteVertices(cgns, m, gvn);
  //
  std::vector<apf::Mesh::Type> apfElementOrder;
  apfElementOrder.push_back(apf::Mesh::QUAD);
  apfElementOrder.push_back(apf::Mesh::TRIANGLE);

  std::vector<CGNS_ENUMT(ElementType_t)> cgnsElementOrder;
  cgnsElementOrder.push_back(CGNS_ENUMV(QUAD_4));
  cgnsElementOrder.push_back(CGNS_ENUMV(TRI_3));
  //
  auto &&cellResult = WriteElements(cgns, m, gvn, 2, faceCount, apfElementOrder, cgnsElementOrder);

  apf::GlobalNumbering *gcn = nullptr;
  gcn = apf::makeGlobal(
      apf::numberOwnedDimension(m, "2D element-nums", 2));
  synchronize(gcn);
  //
  //WriteTags(cgns, std::get<0>(cellResult), std::get<1>(cellResult), std::get<0>(vertResult), std::get<1>(vertResult), std::get<2>(vertResult), m);
  //
  WriteFields(cgns, std::get<0>(cellResult), std::get<1>(cellResult), std::get<0>(vertResult), std::get<1>(vertResult), std::get<2>(vertResult), m);

  destroyGlobalNumbering(gcn);
}

// copy is deliberate
void Write2DEdges(CGNS cgns, apf::Mesh *m, const Count &edgeCount, const Count &vertexCount, apf::GlobalNumbering *gvn)
{
  std::array<cgsize_t, 3> sizes;
  sizes[0] = vertexCount.first; // global
  sizes[1] = edgeCount.first;   // global
  sizes[2] = 0;                 // nodes are unsorted, as defined by api

  {
    std::string baseName("EdgeMeshBase");
    if (cg_base_write(cgns.index, baseName.c_str(), 1, cgns.phys_dim, &cgns.base))
      cg_error_exit();
  }
  // Write the default units at the top.
  if (cg_goto(cgns.index, cgns.base, "end"))
    cg_error_exit();

  if (cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin),
                     CGNS_ENUMV(Degree)))
    cg_error_exit();

  if (cg_dataclass_write(CGNS_ENUMV(Dimensional)))
    cg_error_exit();

  {
    std::string zoneName("EdgeMeshZone");
    if (cg_zone_write(cgns.index, cgns.base, zoneName.c_str(), sizes.data(), CGNS_ENUMV(Unstructured), &cgns.zone))
      cg_error_exit();
  }

  // Stupid api why in the name would you have to repeat ALL vertices again....
  auto &&vertResult = WriteVertices(cgns, m, gvn);
  //
  std::vector<apf::Mesh::Type> apfElementOrder;
  apfElementOrder.push_back(apf::Mesh::EDGE);

  std::vector<CGNS_ENUMT(ElementType_t)> cgnsElementOrder;
  cgnsElementOrder.push_back(CGNS_ENUMV(BAR_2));
  //
  auto &&cellResult = WriteElements(cgns, m, gvn, 1, edgeCount, apfElementOrder, cgnsElementOrder);

  apf::GlobalNumbering *gcn = nullptr;
  gcn = apf::makeGlobal(
      apf::numberOwnedDimension(m, "1D element-nums", 1));
  synchronize(gcn);
  //
  //WriteTags(cgns, std::get<0>(cellResult), std::get<1>(cellResult), std::get<0>(vertResult), std::get<1>(vertResult), std::get<2>(vertResult), m);
  //
  WriteFields(cgns, std::get<0>(cellResult), std::get<1>(cellResult), std::get<0>(vertResult), std::get<1>(vertResult), std::get<2>(vertResult), m);

  destroyGlobalNumbering(gcn);
}

// Todo split this out into a list of calls to local functions to show process/work flow
void WriteCGNS(const char *prefix, apf::Mesh *m, const apf::CGNSBCMap &cgnsBCMap)
{
  static_assert(std::is_same<cgsize_t, int>::value, "cgsize_t not compiled as int");

  const auto myRank = m->getPCU()->Self();
  const Count vertexCount = count(m, 0);
  const Count edgeCount = count(m, 1);
  const Count faceCount = count(m, 2);
  const auto cell_dim = m->getDimension();
  const Count cellCount = count(m, cell_dim);
  // for (int i = 0; i < PCU_Comm_Peers(); ++i)
  // {
  //   if (i == PCU_Comm_Self())
  //   {
  //     std::cout << "*******Local Mesh Stats******************\n";
  //     std::cout << "Rank: " << myRank << ": Number of cells " << cellCount.second << "\n";
  //     std::cout << "Rank: " << myRank << ": Number of vertices " << vertexCount.second << "\n";
  //     std::cout << "Rank: " << myRank << ": Number of faces " << faceCount.second << "\n";
  //     std::cout << "Rank: " << myRank << ": Number of edges " << edgeCount.second << "\n";
  //     std::cout << "*****************************************\n";
  //   }
  //   PCU_Barrier();
  // }

  m->getPCU()->Barrier();
  if (myRank == 0)
  {
    std::cout << "*******Global Mesh Stats*****************\n";
    std::cout << ": Number of cells " << cellCount.first << "\n";
    std::cout << ": Number of vertices " << vertexCount.first << "\n";
    std::cout << ": Number of faces " << faceCount.first << "\n";
    std::cout << ": Number of edges " << edgeCount.first << "\n";
    std::cout << "*****************************************\n";
  }

  std::array<cgsize_t, 3> sizes;
  sizes[0] = vertexCount.first; // global
  sizes[1] = cellCount.first;   // global
  sizes[2] = 0;                 // nodes are unsorted, as defined by api

  // Copy communicator
  auto communicator = m->getPCU()->GetMPIComm();
  cgp_mpi_comm(communicator);
  //
  cgp_pio_mode(CGNS_ENUMV(CGP_INDEPENDENT));

  CGNS cgns;
  cgns.fname = std::string(prefix);
  if (cgp_open(prefix, CGNS_ENUMV(CG_MODE_WRITE), &cgns.index))
    cgp_error_exit();

  {
    std::string baseName("VolumeMeshBase");
    if (cg_base_write(cgns.index, baseName.c_str(), cell_dim, cgns.phys_dim, &cgns.base))
      cg_error_exit();
  }
  // Write the default units at the top.
  if (cg_goto(cgns.index, cgns.base, "end"))
    cg_error_exit();

  if (cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin),
                     CGNS_ENUMV(Degree)))
    cg_error_exit();

  if (cg_dataclass_write(CGNS_ENUMV(Dimensional)))
    cg_error_exit();

  {
    std::string zoneName("VolumeMeshZone");
    if (cg_zone_write(cgns.index, cgns.base, zoneName.c_str(), sizes.data(), CGNS_ENUMV(Unstructured), &cgns.zone))
      cg_error_exit();
  }

  apf::GlobalNumbering *gvn = nullptr;
  gvn = apf::makeGlobal(apf::numberOwnedNodes(m, "node-nums"));
  apf::synchronize(gvn);
  //
  auto &&vertResult = WriteVertices(cgns, m, gvn);
  //
  std::vector<apf::Mesh::Type> apfElementOrder;
  apfElementOrder.push_back(apf::Mesh::HEX);
  apfElementOrder.push_back(apf::Mesh::TET);
  apfElementOrder.push_back(apf::Mesh::PYRAMID);
  apfElementOrder.push_back(apf::Mesh::QUAD);
  apfElementOrder.push_back(apf::Mesh::TRIANGLE);
  apfElementOrder.push_back(apf::Mesh::EDGE);

  std::vector<CGNS_ENUMT(ElementType_t)> cgnsElementOrder;
  cgnsElementOrder.push_back(CGNS_ENUMV(HEXA_8));
  cgnsElementOrder.push_back(CGNS_ENUMV(TETRA_4));
  cgnsElementOrder.push_back(CGNS_ENUMV(PYRA_5));
  cgnsElementOrder.push_back(CGNS_ENUMV(QUAD_4));
  cgnsElementOrder.push_back(CGNS_ENUMV(TRI_3));
  cgnsElementOrder.push_back(CGNS_ENUMV(BAR_2));
  //
  auto &&cellResult = WriteElements(cgns, m, gvn, m->getDimension(), cellCount, apfElementOrder, cgnsElementOrder);
  //
  apf::GlobalNumbering *gcn = nullptr;
  gcn = apf::makeGlobal(apf::numberElements(m, "3D element-nums"));
  // Could be a boco at any dim
  std::map<apf::Mesh::Type, CGNS_ENUMT(ElementType_t)> apf2cgns;
  apf2cgns.insert(std::make_pair(apf::Mesh::HEX, CGNS_ENUMV(HEXA_8)));
  apf2cgns.insert(std::make_pair(apf::Mesh::TET, CGNS_ENUMV(TETRA_4)));
  apf2cgns.insert(std::make_pair(apf::Mesh::PYRAMID, CGNS_ENUMV(PYRA_5)));
  apf2cgns.insert(std::make_pair(apf::Mesh::QUAD, CGNS_ENUMV(QUAD_4)));
  apf2cgns.insert(std::make_pair(apf::Mesh::TRIANGLE, CGNS_ENUMV(TRI_3)));
  apf2cgns.insert(std::make_pair(apf::Mesh::EDGE, CGNS_ENUMV(BAR_2)));
  //
  AddBocosToMainBase(cgns, cellResult, cellCount.first, m, cgnsBCMap, apf2cgns, gvn);
  //
  //WriteTags(cgns, std::get<0>(cellResult), std::get<1>(cellResult), std::get<0>(vertResult), std::get<1>(vertResult), std::get<2>(vertResult), m);
  //
  WriteFields(cgns, std::get<0>(cellResult), std::get<1>(cellResult), std::get<0>(vertResult), std::get<1>(vertResult), std::get<2>(vertResult), m);
  //
  if (cell_dim == 3)
  {
    Write3DFaces(cgns, m, faceCount, vertexCount, gvn);
    //
    Write2DEdges(cgns, m, edgeCount, vertexCount, gvn);
  }
  else if (cell_dim == 2)
  {
    Write2DEdges(cgns, m, edgeCount, vertexCount, gvn);
  }
  //
  destroyGlobalNumbering(gvn);
  destroyGlobalNumbering(gcn);
  //
  cgp_close(cgns.index);
}
} // namespace

namespace apf
{

void writeCGNS(const char *prefix, Mesh *m, const apf::CGNSBCMap &cgnsBCMap)
{
#ifdef HAVE_CGNS
  WriteCGNS(prefix, m, cgnsBCMap);
#else
  PCU_ALWAYS_ASSERT_VERBOSE(true == false,
                            "Build with ENABLE_CGNS to allow this functionality.");
  exit(EXIT_FAILURE);
#endif
}

} // namespace apf
