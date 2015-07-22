#ifndef CHEF_H
#define CHEF_H

#include <apf.h>
#include <apfMesh2.h>
#include <gmi.h>

namespace chef {
  struct IStream;
  struct OStream;
  /** @brief read and write to and from files */
  void cook(gmi_model*& g, apf::Mesh2*& m);
  /** @brief read from stream and write to files */
  void cook(gmi_model*& g, apf::Mesh2*& m, IStream* in);
  /** @brief read from files and write to stream */
  void cook(gmi_model*& g, apf::Mesh2*& m, OStream* out);
  /** @brief read and write to and from streams */
  //void cook(gmi_model*& g, apf::Mesh2*& m, IStream* in, OStream* out);
  /** @brief create an input stream from an output stream */
  IStream* makeIStream(OStream* ostream);
  /** @brief destroy input stream */
  void destroyIStream(IStream* is);
  /** @brief make output stream */
  OStream* makeOStream();
  /** @brief destroy output stream */
  void destroyOStream(OStream* os);
}

#endif
