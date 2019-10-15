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

#ifdef HAVE_CGNS
//
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include <utility>
//
#include <cgns_io.h>
#include <pcgnslib.h>
//
#endif

namespace
{
#ifdef DEBUG
static constexpr bool debugOutput = true;
#else
static constexpr bool debugOutput = false;
#endif

static std::string CGNSElementTypeToString(const CGNS_ENUMT(ElementType_t) & elementType)
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
  if constexpr (debugOutput)
  {
    for (int i = 0; i < PCU_Comm_Peers(); i++)
    {
      if (i == PCU_Comm_Self())
      {
        out << "Rank [" << i << "] " << std::forward<Arg>(arg);
        ((out << ", " << std::forward<Args>(args)), ...);
        out << "\n";
        out << std::flush;
      }
      PCU_Barrier();
    }
  }
}

void Kill()
{
  if (PCU_Comm_Initialized())
  {
    // Finalize the MPI environment.
    PCU_Comm_Free();
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  else
  {
    exit(EXIT_FAILURE);
  }
}

auto ReadCGNSCoords(int cgid, int base, int zone, int ncoords, int nverts, const std::vector<cgsize_t> &, const apf::GlobalToVert &globalToVert)
{
  // Read all
  // const int lowest = 1;
  // const int highest = nverts;

  // Read min required as defined by consecutive range
  // make one based as ReadElements makes zero based
  // const int lowest = *std::min_element(vertexIDs.begin(), vertexIDs.end()) + 1;
  // const int highest = *std::max_element(vertexIDs.begin(), vertexIDs.end()) + 1;
  // DebugParallelPrinter(std::cout, "From vertexIDs ", lowest, highest, nverts);

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
      Kill();
    }
    const auto coord_name = std::string(coord_names[d].c_str());
    //boost::algorithm::trim(coord_name); // can't be bothered including boost

    auto &ordinate = ordinates[d];
    if (cgp_coord_read_data(cgid, base, zone, d + 1, &range_min[0], &range_max[0],
                            ordinate.data()))
    {
      std::cout << __LINE__ << " CGNS is dead " << std::endl;
      Kill();
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

apf::Mesh2 *DoIt(gmi_model* g, const std::string &fname, int readDim = 0)
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
    Kill();
  }

  std::string basename;
  basename.resize(CGIO_MAX_NAME_LENGTH + 1, ' ');

  int cellDim = -1;
  int physDim = -1;
  const int base = 1;
  cg_base_read(cgid, base, &basename[0], &cellDim, &physDim);
  if (readDim == 0)
    readDim = cellDim;

  // Salome cgns is a bit on the piss: cellDim, physDim, ncoords are not always consistent
  apf::Mesh2 *mesh = apf::makeEmptyMdsMesh(g, cellDim, false);
  apf::GlobalToVert globalToVert;

  int nzones = -1;
  cg_nzones(cgid, base, &nzones);
  basename = std::string(basename.c_str());

  int nfam = -1;
  cg_nfamilies(cgid, base, &nfam);

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
      Kill();
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
      Kill();
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

      CGNSElementTypeToString(elementType);
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
            mesh->setPoint(globalToVert[p.first], 0, point);
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
                  << " " << CGNSElementTypeToString(elementType) << std::endl;
        Kill();
      }
    }
  }
  cgp_close(cgid);

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

  return mesh;
}
} // namespace

namespace apf
{

// caller needs to bring up and pull down mpi/pcu: mpi/pcu is required and assumed.
Mesh2 *loadMdsFromCGNS(gmi_model* g, const char *fname)
{
#ifdef HAVE_CGNS
  Mesh2 *m = DoIt(g, fname);
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
