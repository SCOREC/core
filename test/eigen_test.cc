#include <apfMatrix.h>
#include <algorithm>

struct Input {
  double A[3][3];
  double V[3][3];
  double l[3];
};

#define NINPUTS 6

static struct Input const inputs[NINPUTS] = {
{{{1.001575e+00, -3.138397e-01, 8.107355e-01},
  {-3.138397e-01, 4.946182e-01, -1.860431e+00},
  {8.107355e-01, -1.860431e+00, 7.582283e+00}},
 {{1.238194e-01, -9.850387e-01, 1.198649e-01},
  {9.665579e-01, 9.238558e-02, -2.392296e-01},
  {2.245766e-01, 1.454777e-01, 9.635360e-01}},
 {2.214902e-02, 9.112746e-01, 8.145053e+00}},
{{{1.668786e+00, -7.051850e-02, 4.561114e-01},
  {-7.051850e-02, 1.022233e+00, -5.465726e-01},
  {4.561114e-01, -5.465726e-01, 1.128792e+00}},
 {{-2.295847e-01, 6.324411e-01, -7.398034e-01},
  {6.673195e-01, 6.556027e-01, 3.533692e-01},
  {7.085023e-01, -4.125571e-01, -5.725566e-01}},
 {4.661900e-01, 1.298152e+00, 2.055468e+00}},
{{{7.139854e-01, 1.220499e+00, -8.923666e-02},
  {1.220499e+00, 2.289174e+00, -4.099688e-01},
  {-8.923666e-02, -4.099688e-01, 6.159893e-01}},
 {{8.433346e-01, 2.675115e-01, -4.660733e-01},
  {-4.929524e-01, 3.976489e-02, -8.691471e-01},
  {-2.139735e-01, 9.627338e-01, 1.654055e-01}},
 {2.321140e-02, 5.742601e-01, 3.021677e+00}},
{{{1.662534e+00, 1.834307e+00, 1.113124e-01},
  {1.834307e+00, 4.465579e+00, 2.825547e+00},
  {1.113124e-01, 2.825547e+00, 3.277306e+00}},
 {{-6.417678e-01, -7.179301e-01, 2.696488e-01},
  {5.861286e-01, -2.324275e-01, 7.761641e-01},
  {-4.945578e-01, 6.561660e-01, 5.699638e-01}},
 {7.303538e-02, 2.154649e+00, 7.177734e+00}},
{{{2.869781e+00, -1.879048e+00, 1.982487e+00},
  {-1.879048e+00, 1.296580e+00, -9.613331e-01},
  {1.982487e+00, -9.613331e-01, 5.662202e+00}},
 {{5.819343e-01, -6.425814e-01, 4.984392e-01},
  {8.105165e-01, 5.083603e-01, -2.909171e-01},
  {-6.644880e-02, 5.732879e-01, 8.166551e-01}},
 {2.627564e-02, 2.587634e+00, 7.214654e+00}},
{{{2.763219e+00, 1.444043e+00, 1.102805e+00},
  {1.444043e+00, 4.883458e+00, 4.285960e+00},
  {1.102805e+00, 4.285960e+00, 5.790849e+00}},
 {{-2.258181e-01, 9.448099e-01, 2.373612e-01},
  {7.543796e-01, 1.543447e-02, 6.562570e-01},
  {-6.163746e-01, -3.272552e-01, 7.162307e-01}},
 {9.493002e-01, 2.404829e+00, 1.008340e+01}}
};

static void toVectors(apf::Matrix3x3 ev, apf::Vector3 ew,
    apf::Vector<4> out[3])
{
  for (int i = 0; i < 3; ++i) {
    out[i][0] = ew[i];
    for (int j = 0; j < 3; ++j)
      out[i][j + 1] = ev[i][j];
  }
}

static void fromVectors(apf::Matrix3x3 & ev, apf::Vector3 & ew,
    apf::Vector<4> out[3])
{
  for (int i = 0; i < 3; ++i) {
    ew[i] = out[i][0];
    for (int j = 0; j < 3; ++j)
      ev[i][j] = out[i][j + 1];
  }
}

struct VectorLess {
  bool operator()(apf::Vector<4> const& a, apf::Vector<4> const& b)
  {
    return a[0] < b[0];
  }
};

static void sortVectors(apf::Vector<4> v[3])
{
  std::sort(v, v + 3, VectorLess());
}

static void sortEigen(apf::Matrix3x3 & ev, apf::Vector3 & ew)
{
  apf::Vector<4> v[3];
  toVectors(ev, ew, v);
  sortVectors(v);
  fromVectors(ev, ew, v);
}

static double diffnorm(double a, double b)
{
  double diff = a - b;
  return diff * diff;
}

static double diffnorm(apf::Vector<3> const& a, apf::Vector<3> const& b)
{
  double n = 0;
  for (int i = 0; i < 3; ++i)
    n += diffnorm(a[i], b[i]);
  return n;
}

static double diffnormEigenVec(apf::Vector<3> const& a, apf::Vector<3> const& b)
{
  /* may be pointing in the opposite direction, correct for that */
  if (a * b < 0)
    return diffnorm(a, b * -1);
  return diffnorm(a, b);
}

static double diffnormEigenVecs(apf::Matrix3x3 const& a, apf::Matrix3x3 const& b)
{
  double n = 0;
  for (int i = 0; i < 3; ++i)
    n += diffnormEigenVec(a[i], b[i]);
  return n;
}

int main()
{
  for (int i = 0; i < NINPUTS; ++i) {
    apf::Matrix3x3 A(&(inputs[i].A[0]));
    apf::Matrix3x3 V(&(inputs[i].V[0]));
    V = apf::transpose(V); /* Octave generated eigenvector columns, we generate rows */
    apf::Vector3 l(&(inputs[i].l[0]));
    apf::Matrix3x3 V2;
    apf::Vector3 l2;
    int n = apf::eigen(A, &V2[0], &l2[0]);
    assert(n == 3);
    sortEigen(V, l);
    sortEigen(V2, l2);
    assert(diffnormEigenVecs(V, V2) < 1e-10);
    assert(diffnorm(l, l2) < 1e-10);
  }
}
