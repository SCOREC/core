#ifndef APFCOMPLEX_H_
#define APFCOMPLEX_H_

#ifdef C_COMPLEX
  #include <complex.h>
#endif

#define CXX_COMPLEX 1
#ifdef CXX_COMPLEX
  #include <complex>
  using double_complex = std::complex<double>;
#endif

#endif
