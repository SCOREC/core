/*
 * Author: Dr Andrew Parker (2019) - FGE Ltd
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include "apfMesh.h"
#include "apfNumbering.h"
#include "apfNumberingClass.h"
#include "apfShape.h"
#include "apfFieldData.h"
#include <pcu_util.h>
#include <lionPrint.h>
//
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

namespace
{
//https://proteustoolkit.org/capi/html/_parallel_mesh_converter_8cpp_source.html
static auto count(apf::Mesh *m, int dim)
{
  const int local = apf::countOwned(m, dim);
  int total = local;
  PCU_Add_Ints(&total, 1); // size of total array
  return std::make_pair(total, local);
}

void WriteCGNS(const char *prefix, apf::Mesh *m, const apf::CGNSBCMap &cgnsBCMap)
{
  //ShowNumbering(m);

  const auto myRank = PCU_Comm_Self();
  const auto vertexCount = count(m, 0);
  const auto edgeCount = count(m, 1);
  const auto faceCount = count(m, 2);
  const auto cell_dim = m->getDimension();
  const auto cellCount = count(m, cell_dim);
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

  PCU_Barrier();
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
  sizes[0] = vertexCount.first;
  sizes[1] = cellCount.first;
  sizes[2] = 0; // nodes are unsorted.

  // Copy communicator
  auto communicator = PCU_Get_Comm();
  cgp_mpi_comm(communicator);
  //
  cgp_pio_mode(CGNS_ENUMV(CGP_INDEPENDENT));

  int index = -1;
  if (cgp_open(prefix, CGNS_ENUMV(CG_MODE_WRITE), &index))
    cgp_error_exit();

  int base = -1;
  const int phys_dim = 3;
  {
    std::string baseName("Base_" + std::to_string(1));
    if (cg_base_write(index, baseName.c_str(), cell_dim, phys_dim, &base))
      cg_error_exit();
  }
  // Write the default units at the top.
  if (cg_goto(index, base, "end"))
    cg_error_exit();

  if (cg_units_write(CGNS_ENUMV(Kilogram), CGNS_ENUMV(Meter), CGNS_ENUMV(Second), CGNS_ENUMV(Kelvin),
                     CGNS_ENUMV(Degree)))
    cg_error_exit();

  if (cg_dataclass_write(CGNS_ENUMV(Dimensional)))
    cg_error_exit();

  int zone = -1;
  {
    std::string zoneName("Zone_" + std::to_string(1));
    if (cg_zone_write(index, base, zoneName.c_str(), sizes.data(), CGNS_ENUMV(Unstructured), &zone))
      cg_error_exit();
  }

  static_assert(std::is_same<cgsize_t, int>::value, "cgsize_t not compiled as int");

  int Cx = -1;
  int Cy = -1;
  int Cz = -1;

  if (phys_dim > 0)
  {
    if (cgp_coord_write(index, base, zone, CGNS_ENUMV(RealDouble), "CoordinateX", &Cx))
      cgp_error_exit();
  }
  if (phys_dim > 1)
  {
    if (cgp_coord_write(index, base, zone, CGNS_ENUMV(RealDouble), "CoordinateY", &Cy))
      cgp_error_exit();
  }
  if (phys_dim > 2)
  {
    if (cgp_coord_write(index, base, zone, CGNS_ENUMV(RealDouble), "CoordinateZ", &Cz))
      cgp_error_exit();
  }

  std::array<std::vector<double>, 3> coords;

  apf::GlobalNumbering *gvn = nullptr;
  gvn = apf::makeGlobal(apf::numberOwnedNodes(m, "node-nums"));
  apf::synchronize(gvn);

  cgsize_t rmin[3];
  cgsize_t rmax[3];
  rmin[0] = std::numeric_limits<cgsize_t>::max();
  rmax[0] = 0;

  {
    apf::Vector3 point;
    for (int i = 0; i < PCU_Comm_Peers(); ++i)
    {
      if (i == PCU_Comm_Self())
      {
        apf::MeshIterator *vertIter = m->begin(0);
        apf::MeshEntity *vert = nullptr;
        while ((vert = m->iterate(vertIter)))
        {
          if (m->isOwned(vert))
          {
            const cgsize_t n = static_cast<cgsize_t>(apf::getNumber(gvn, vert, 0) + 1); // one based
            rmin[0] = std::min(rmin[0], n);
            rmax[0] = std::max(rmax[0], n);

            m->getPoint(vert, 0, point);
            coords[0].push_back(point[0]);
            coords[1].push_back(point[1]);
            coords[2].push_back(point[2]);
          }
        }
        m->end(vertIter);
      }
    }
  }

  // oddness of the api
  rmin[1] = rmin[0];
  rmin[2] = rmin[0];
  rmax[1] = rmax[0];
  rmax[2] = rmax[0];

  if (phys_dim > 0)
  {
    if (cgp_coord_write_data(index, base, zone, Cx, &rmin[0], &rmax[0], coords[0].data()))
      cgp_error_exit();
  }
  if (phys_dim > 1)
  {
    if (cgp_coord_write_data(index, base, zone, Cy, &rmin[0], &rmax[0], coords[1].data()))
      cgp_error_exit();
  }
  if (phys_dim > 2)
  {
    if (cgp_coord_write_data(index, base, zone, Cz, &rmin[0], &rmax[0], coords[2].data()))
      cgp_error_exit();
  }

  apf::GlobalNumbering *gcn = nullptr;
  gcn = apf::makeGlobal(apf::numberElements(m, "element-nums"));

  std::vector<int> apfElementOrder;
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

  std::vector<int> globalNumbersByElementType(apfElementOrder.size(), 0);
  std::vector<int> numbersByElementType(apfElementOrder.size(), 0);
  for (std::size_t o = 0; o < apfElementOrder.size(); o++)
  {
    apf::MeshIterator *cellIter = m->begin(cell_dim);
    apf::MeshEntity *cell = nullptr;
    int counter = 0;
    while ((cell = m->iterate(cellIter)))
    {
      if (m->getType(cell) == apfElementOrder[o] && m->isOwned(cell))
      {
        counter++;
      }
    }
    m->end(cellIter);
    numbersByElementType[o] = counter;
    int total = counter;
    PCU_Add_Ints(&total, 1); // size of total array
    globalNumbersByElementType[o] = total;
  }
  cgsize_t allTotal = std::accumulate(globalNumbersByElementType.begin(), globalNumbersByElementType.end(), 0);
  PCU_ALWAYS_ASSERT_VERBOSE(allTotal == cellCount.first, ("Must be equal " + std::to_string(allTotal) + " " + std::to_string(cellCount.first)).c_str());

  int globalStart = 1; // one-based
  for (std::size_t o = 0; o < apfElementOrder.size(); o++)
  {
    std::vector<cgsize_t> elements;
    apf::MeshIterator *cellIter = m->begin(cell_dim);
    apf::MeshEntity *cell = nullptr;
    apf::Downward verts;

    while ((cell = m->iterate(cellIter)))
    {
      if (m->getType(cell) == apfElementOrder[o] && m->isOwned(cell))
      {
        const auto numVerts = m->getDownward(cell, 0, verts);
        for (int i = 0; i < numVerts; i++)
        {
          const auto n = apf::getNumber(gvn, verts[i], 0);
          elements.push_back(n + 1); // one-based
        }
      }
    }
    m->end(cellIter);

    if (globalNumbersByElementType[o] > 0)
    {
      const int globalEnd = globalStart + globalNumbersByElementType[o] - 1; // one-based stuff
      //
      int sectionNumber = -1;
      if (cgp_section_write(index, base, zone, (std::string(cg_ElementTypeName(cgnsElementOrder[o])) + " " + std::to_string(globalStart) + "->" + std::to_string(globalEnd)).c_str(), cgnsElementOrder[o], globalStart,
                            globalEnd, 0, &sectionNumber)) // global start, end within the file for that element type
        cgp_error_exit();

      std::vector<int> allNumbersForThisType(PCU_Comm_Peers(), 0);
      MPI_Allgather(&numbersByElementType[o], 1, MPI_INT, allNumbersForThisType.data(), 1,
                    MPI_INT, PCU_Get_Comm());

      cgsize_t num = 0;
      for (int i = 0; i < PCU_Comm_Self(); i++)
        num += allNumbersForThisType[i];

      cgsize_t elStart = globalStart + num;
      cgsize_t elEnd = elStart + numbersByElementType[o] - 1;                       // one-based stuff
      if (cgp_elements_write_data(index, base, zone, sectionNumber, elStart, elEnd, // per processor within the range[start, end]
                                  elements.data()))
        cgp_error_exit();

      //std::cout << "RANK: " << PCU_Comm_Self() << " ==> " << globalStart << " " << globalEnd << " elStart " << elStart << " elEnd " << elEnd << " numbersByElementType[o] " << numbersByElementType[o] << std::endl;

      globalStart += globalNumbersByElementType[o];
    }
  }
  //
  std::cout << &cgnsBCMap << std::endl;
  if (cell_dim == 3)
  {
    auto iter = cgnsBCMap.find("Vertex");
    if (iter != cgnsBCMap.end())
    {
      for (const auto &p : iter->second)
      {
        std::cout << iter->first << " " << p.first << " " << p.second << std::endl;
      }
    }
    iter = cgnsBCMap.find("EdgeCenter");
    if (iter != cgnsBCMap.end())
    {
      for (const auto &p : iter->second)
      {
        std::cout << iter->first << " " << p.first << " " << p.second << std::endl;
      }
    }
    iter = cgnsBCMap.find("FaceCenter");
    if (iter != cgnsBCMap.end())
    {
      for (const auto &p : iter->second)
      {
        std::cout << iter->first << " " << p.first << " " << p.second << std::endl;
        auto* field = m->findField(p.first.c_str());
        std::cout << field << " " << p.second << std::endl;

      }
    }
    iter = cgnsBCMap.find("CellCenter");
    if (iter != cgnsBCMap.end())
    {
      for (const auto &p : iter->second)
      {
        std::cout << iter->first << " " << p.first << " " << p.second << std::endl;
      }
    }
  }
  //
  destroyGlobalNumbering(gvn);
  destroyGlobalNumbering(gcn);
  //
  cgp_close(index);
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
