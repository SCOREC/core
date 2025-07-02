#include "apfCAP.h"

#include <exception>
#include <fstream>
#include <stdexcept>

#include <apf.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <pcu_util.h>

#include <CreateMG_Framework_Mesh.h>
#ifdef PUMI_HAS_CAPSTONE_SIZINGMETRICTOOL
#include <CreateMG_SizingMetricTool.h>
#endif

namespace {

/** \brief Load bulk file of one type into std::vector. */
template <class T>
typename std::enable_if<std::is_trivially_copyable<T>::value, void>::type
load_file(const std::string& fname, std::vector<T>& data) {
  try {
    std::ifstream file(fname, std::ios::in | std::ios::binary);
    if (!file.is_open()) {
      throw std::runtime_error("failed to open file `" + fname + "`");
    }
    T val;
    while (file) {
      file.read(reinterpret_cast<char*>(std::addressof(val)), sizeof(T));
      if (file.gcount() == sizeof(T)) data.push_back(val);
    }
  } catch (...) {
    std::throw_with_nested(
      std::runtime_error("failed to load file `" + fname + "`")
    );
  }
}

bool doSmoothing(
  apf::Mesh2* m, const std::string& analysis,
  const std::vector<CreateMG::Metric6>& in, std::vector<CreateMG::Metric6>& out
) {
  auto smooth_tool = CreateMG::get_sizing_metric_tool(
    apf::exportCapNative(m)->get_context(), "CreateSmoothingBase"
  );
  if (smooth_tool == nullptr) {
    lion_eprint(1, "ERROR: Unable to find \"CreateSmoothingBase\"\n");
    return false;
  }
  smooth_tool->set_context(apf::exportCapNative(m)->get_context());
  CreateMG::M_MModel mmodel;
  MG_API_CALL(apf::exportCapNative(m), get_current_model(mmodel));
  out.resize(in.size());
  smooth_tool->set_metric(mmodel, "apf_cap_sizing6", in);
  smooth_tool->smooth_metric(mmodel, analysis, "apf_cap_sizing6", out);
  smooth_tool->remove(mmodel, "apf_cap_sizing6");
  return true;
}

} // namespace

namespace apf {

bool loadCapSizing(
  apf::Mesh2* m, const std::vector<CreateMG::Metric6>& sizing,
  apf::Field* scales, apf::Field* frames
) {
  PCU_DEBUG_ASSERT(m);
  PCU_DEBUG_ASSERT(scales);
  PCU_DEBUG_ASSERT(frames);
  apf::Vector3   scale;
  apf::Matrix3x3 frame;
  for (std::size_t i = 0; i < sizing.size(); ++i) {
    // field eigen values
    apf::Matrix3x3 mat(
      sizing[i][0], sizing[i][1], sizing[i][2],
      sizing[i][1], sizing[i][3], sizing[i][4],
      sizing[i][2], sizing[i][4], sizing[i][5]
    );
    int n = apf::eigen(mat, &frame[0], &scale[0]);
    PCU_DEBUG_ASSERT(n == 3);
    for (int i = 0; i < 3; i++) {
      scale[i] = 1.0/sqrt(scale[i]);
    }
    // MeshAdapt wants frames on the columns.
    frame = apf::transpose(frame);
    apf::MeshEntity* e = getCapEntity(m, 0, i + 1);
    apf::setVector(scales, e, 0, scale);
    apf::setMatrix(frames, e, 0, frame);
  }
  return true;
}

bool loadCapSizing(
  apf::Mesh2* m, const std::vector<CreateMG::Metric6>& sizing,
  const char* scales, const char* frames
) {
  apf::Field* scalesField = apf::createFieldOn(m, scales, VECTOR);
  apf::Field* framesField = apf::createFieldOn(m, frames, MATRIX);
  return loadCapSizing(m, sizing, scalesField, framesField);
}

void loadCapSizingFileMetrics(
  apf::Mesh2* m, const std::string& sizingFile, const std::string& vmapFile,
  std::vector<CreateMG::Metric6>& sizing6
) {
  std::vector<size_t> vmap;
  std::vector<double> sizing;
  load_file(vmapFile, vmap);
  load_file(sizingFile, sizing);
  // Validate sizing
  if (sizing.size() / 6 != m->count(0)) {
    std::string msg = "loadCapSizingFile: sizing file size does not match mesh"
      " vertex count: " + std::to_string(sizing.size() / 6) + " / " +
      std::to_string(m->count(0));
    fail(msg.c_str());
  }
  // Validate vmap
  if (
    vmap.size() < 1 || vmap[0] != vmap.size() - 1 || vmap[0] != m->count(0)
  ) {
    fail("loadCapSizingFile: loaded invalid vmap");
  }
  // Copy sizefield to Metric6
  sizing6.resize(vmap[0]);
  for (size_t si = 0, vi = 1; vi < vmap.size(); ++vi) {
    size_t tid = vmap[vi];
    PCU_DEBUG_ASSERT(tid != 0);
    for (size_t j = 0; j < 6; ++j, ++si) {
      sizing6[tid - 1][j] = sizing[si];
    }
  }
}

bool loadCapSizingFile(
  apf::Mesh2* m, const std::string& sizingFile, const std::string& vmapFile,
  apf::Field* scales, apf::Field* frames,
  bool smooth, const std::string& analysis
) {
  if (smooth && !has_smoothCapAnisoSizes()) {
    lion_eprint(1,
      "WARNING: loadCapSizingFile: smoothing requested but apf_cap was"
      " compiled without support for smoothing\n");
    return false;
  }
  std::vector<CreateMG::Metric6> sizing6;
  loadCapSizingFileMetrics(m, sizingFile, vmapFile, sizing6);
  if (smooth) {
    std::vector<CreateMG::Metric6> ometric;
    if (!doSmoothing(m, analysis, sizing6, ometric)) return false;
    std::swap(sizing6, ometric);
  }
  return loadCapSizing(m, sizing6, scales, frames);
}

bool loadCapSizingFile(
  apf::Mesh2* m, const std::string& sizingFile, const std::string& vmapFile,
  const char* scales, const char* frames,
  bool smooth, const std::string& analysis
) {
  apf::Field* scalesField = apf::createFieldOn(m, scales, VECTOR);
  apf::Field* framesField = apf::createFieldOn(m, frames, MATRIX);
  return loadCapSizingFile(
    m, sizingFile, vmapFile, scalesField, framesField, smooth, analysis
  );
}

void extractCapSizing(
  apf::Mesh2* mesh, apf::Field* scales, apf::Field* frames,
  std::vector<CreateMG::Metric6>& sizing
) {
  apf::Matrix3x3 Q;
  apf::Vector3 H;
  sizing.resize(mesh->count(0));
  apf::MeshIterator* it = mesh->begin(0);
  for (apf::MeshEntity* e = mesh->iterate(it); e; e = mesh->iterate(it)) {
    apf::getVector(scales, e, 0, H); // Desired element lengths.
    apf::getMatrix(frames, e, 0, Q); // MeshAdapt uses column vectors.
    apf::Matrix3x3 L(
      1.0/(H[0]*H[0]), 0, 0,
      0, 1.0/(H[1]*H[1]), 0,
      0, 0, 1.0/(H[2]*H[2])
    );
    apf::Matrix3x3 t = Q * L * apf::transpose(Q); // Invert orthogonal frames.
    size_t id = getCapId(mesh, e);
    PCU_DEBUG_ASSERT(id != 0);
    --id;
    sizing[id][0] = t[0][0];
    sizing[id][1] = t[0][1];
    sizing[id][2] = t[0][2];
    sizing[id][3] = t[1][1];
    sizing[id][4] = t[1][2];
    sizing[id][5] = t[2][2];
  }
  mesh->end(it);
}

bool has_smoothCapAnisoSizes(void) noexcept {
#ifdef PUMI_HAS_CAPSTONE_SIZINGMETRICTOOL
  return true;
#else
  return false;
#endif
}

bool smoothCapAnisoSizes(apf::Mesh2* mesh, std::string analysis,
  apf::Field* scales, apf::Field* frames) {
#ifdef PUMI_HAS_CAPSTONE_SIZINGMETRICTOOL
  std::vector<CreateMG::Metric6> sizing6, ometric;
  extractCapSizing(mesh, scales, frames, sizing6);
  if (!doSmoothing(mesh, analysis, sizing6, ometric)) return false;
  return loadCapSizing(mesh, ometric, frames, scales);
#else
  (void) mesh;
  (void) analysis;
  (void) scales;
  (void) frames;
  apf::fail("smoothCAPAnisoSizes: Capstone does not have SizingMetricTool.");
#endif
}

} // namespace apf
