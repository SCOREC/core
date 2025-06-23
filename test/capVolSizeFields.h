#include <lionPrint.h>
#include <ma.h>

// analytic size field classes //
class UniformAniso : public ma::AnisotropicFunction
{
  public:
    UniformAniso(ma::Mesh* m)
    {
      mesh = m;
    }
    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    {
      (void)v;
      ma::Vector h;
      // principal scales
      h = ma::Vector(1,
		     0.1,
		     0.1);
      // principal directions
      ma::Matrix r(1.0, 0.0, 0.0,
		   0.0, 1.0, 0.0,
		   0.0, 0.0, 1.0);
      H = h;
      R = r;
    }
  private:
    ma::Mesh* mesh;
};

class WingShock : public ma::AnisotropicFunction
{
  public:
    WingShock(ma::Mesh* m, double inFactor)
    {
      mesh = m;
      factor = inFactor;
    }
    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    {
      ma::Vector p = ma::getPosition(mesh,v);
      double x = p[0];
      double x0 = 0.5;
      double planeSize = 0.03125;
      double spanSize  = 0.5;
      double delta = 0.5;

      double beta = 0.3;
      double x1   = 0.2;
      double f    = beta + x * (1. - beta) / x1;
      double multipier = (x <= x1) ? 1.0 : f;

      double s0 = planeSize / factor;
      double alpha = planeSize * (1. - 1./factor) / delta;
      double hx = s0 + alpha * std::abs(x - x0);
      double hy = spanSize;
      double hz = spanSize;
      hx *= multipier;
      hy *= multipier;

      ma::Vector h;

      h = ma::Vector(hx, hy, hz);
      /* // principal directions */
      ma::Matrix r(1.0, 0.0, 0.0,
		   0.0, 1.0, 0.0,
		   0.0, 0.0, 1.0);
      H = h;
      R = r;
    }
  private:
    ma::Mesh* mesh;
    double factor;
};

class Shock : public ma::AnisotropicFunction
{
  public:
    Shock(ma::Mesh* m)
    {
      mesh = m;
    }
    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    {
      ma::Vector p = ma::getPosition(mesh,v);
      double x = p[0];
      double x0 = 0.0;
      double s = 0.25 * (1 - exp(-2*std::abs(x - x0))) + 0.001;

      double hx = s;
      double hy = 0.25;
      double hz = 0.25;

      ma::Vector h;

      h = ma::Vector(hx, hy, hz);
      /* // principal directions */
      ma::Matrix r(1.0, 0.0, 0.0,
		   0.0, 1.0, 0.0,
		   0.0, 0.0, 1.0);
      H = h;
      R = r;
    }
  private:
    ma::Mesh* mesh;
};

// Shock along a cylinder like the one in
// https://www.scorec.rpi.edu/~cwsmith/SCOREC//SimModSuite_2023.1-230428dev/MeshSimAdapt/group__MSAEx__GenAnisoMesh.html
class CylBoundaryLayer : public ma::AnisotropicFunction {
  ma::Mesh* m;
  enum Axis : int {X = 0, Y = 1, Z = 2} axis;
  double radius, length;

public:
  CylBoundaryLayer(ma::Mesh* mesh): m(mesh) {
    ma::Vector bmin(HUGE_VAL, HUGE_VAL, HUGE_VAL),
               bmax(-HUGE_VAL, -HUGE_VAL, -HUGE_VAL);
    ma::getBoundingBox(mesh, bmin, bmax);
    ma::Vector len = bmax - bmin;
    if (len.x() > len.y() && len.x() > len.z()) {
      lion_oprint(1, "Shock2 along X axis.\n");
      length = len.x();
      axis = Axis::X;
      // Radial lengths should be the same, but we find average diameter and
      // then halve that.
      radius = (len.y() / 2 + len.z() / 2) / 2;
    } else if (len.y() > len.x() && len.y() > len.z()) {
      lion_oprint(1, "Shock2 along Y axis.\n");
      length = len.y();
      axis = Axis::Y;
      radius = (len.x() / 2 + len.z() / 2) / 2;
    } else {
      lion_oprint(1, "Shock2 along Z axis.\n");
      length = len.z();
      axis = Axis::Z;
      radius = (len.x() / 2 + len.y() / 2) / 2;
    }
  }

  void getValue(ma::Entity *v, ma::Matrix &R, ma::Vector &H) {
    const double M = 2.5, aspect_ratio = 16;
    ma::Vector p = ma::getPosition(m, v);
    H = ma::Vector(1, 1, 1) * radius / M;
    R[axis] = ma::Vector(0, 0, 0);
    R[axis][axis] = 1;

    int cosDirs[] = {1, 2, 0}, sinDirs[] = {2, 0, 1}; // Rotated sin/cos axes
    int cosDir = cosDirs[axis], sinDir = sinDirs[axis];

    double distFromAxis = std::sqrt(p[cosDir] * p[cosDir] + p[sinDir] * p[sinDir]);

    if (std::fabs(radius * radius - p[cosDir] * p[cosDir] -
                  p[sinDir] * p[sinDir]) < radius*radius / 4) {
      H[cosDir] /= aspect_ratio;
    }

    double cosT = p[cosDir] / distFromAxis;
    double sinT = p[sinDir] / distFromAxis;
    R[cosDir][cosDir] = cosT; R[cosDir][sinDir] = -sinT; R[cosDir][axis] = 0;
    R[sinDir][cosDir] = sinT; R[sinDir][sinDir] = cosT; R[sinDir][axis] = 0;

    // Normalize R rows. Skip axis which is already normal.
    for (int i = 0; i < 3; ++i) {
      if (i != axis) R[i] = R[i].normalize();
    }
  }
};

class Uniform : public ma::IsotropicFunction
{
  public:
    Uniform(ma::Mesh* m, double inFactor)
    {
      mesh = m;
      factor = inFactor;
      average = ma::getAverageEdgeLength(m);
      printf("average is %f\n", average);
    }
    virtual double getValue(ma::Entity* v)
    {
      (void)v;
      return 0.03125;
    }
  private:
    ma::Mesh* mesh;
    double factor;
    double average;
};

class Linear1 : public ma::IsotropicFunction
{
  public:
    Linear1(ma::Mesh* m, double inFactor)
    {
      mesh = m;
      factor = inFactor;
      average = ma::getAverageEdgeLength(m);
      printf("average is %f\n", average);
    }
    virtual double getValue(ma::Entity* v)
    {
      ma::Vector p = ma::getPosition(mesh, v);
      double x = p[0];
      if (x <= 0)
      	return 5.0;
      else if (x > 0 && x <= 1000)
      	return 10.0;
      else if (x > 1000 && x <= 2000)
      	return 20.0;
      else
      	return 30.0;
    }
  private:
    ma::Mesh* mesh;
    double factor;
    double average;
};

class Linear2 : public ma::IsotropicFunction
{
  public:
    Linear2(ma::Mesh* m, double inFactor)
    {
      mesh = m;
      factor = inFactor;
      average = ma::getAverageEdgeLength(m);
      printf("average is %f\n", average);
    }
    virtual double getValue(ma::Entity* v)
    {
      ma::Vector p = ma::getPosition(mesh, v);
      double x = p[0];
      if (x >= 0 && x <= 20)
      	return 0.5;
      else if (x > 20 && x <= 60)
      	return 2.0;
      else if (x > 60 && x <= 100)
      	return 3.0;
      else
      	return 4.0;
    }
  private:
    ma::Mesh* mesh;
    double factor;
    double average;
};

class Linear3 : public ma::IsotropicFunction
{
  public:
    Linear3(ma::Mesh* m, double inFactor)
    {
      mesh = m;
      factor = inFactor;
      average = ma::getAverageEdgeLength(m);
      ma::getBoundingBox(m, lower, upper);
      printf("average is %f\n", average);
    }
    virtual double getValue(ma::Entity* v)
    {
      ma::Vector p = ma::getPosition(mesh, v);
      double x = p[0];

      // get current size
      apf::Up up;
      mesh->getUp(v, up);
      double sum = 0;
      double count = 0;
      for (int i = 0; i < up.n; i++) {
	ma::Entity* currentEdge = up.e[i];
	ma::Entity* vs[2];
	mesh->getDownward(currentEdge, 0, vs);
	ma::Vector diff = ma::getPosition(mesh, vs[0]) - ma::getPosition(mesh, vs[1]);
	sum += diff.getLength();
	count++;
      }
      double currentSize = sum/count;

      // compute the size multiplier
      double a = (factor/2 - 2/factor) / (upper[0] - lower[0]) / (upper[0] - lower[0]);
      double b = 2/factor;
      double y = a * (x - lower[0]) * (x - lower[0]) + b;

      // new size is multiplier * currentSize
      return y*currentSize;
    }
  private:
    ma::Mesh* mesh;
    double factor;
    double average;
    ma::Vector upper;
    ma::Vector lower;
};

class GeomB737 : public ma::IsotropicFunction
{
  public:
    GeomB737(ma::Mesh* m, double inFactor)
    {
      mesh = m;
      factor = inFactor;
      average = ma::getAverageEdgeLength(m);
      ma::getBoundingBox(m, lower, upper);
      printf("average is %f\n", average);
    }
    virtual double getValue(ma::Entity* v)
    {
      ma::Vector p = ma::getPosition(mesh, v);
      double x = p[0];
      (void)x;

      // get current size
      apf::Up up;
      mesh->getUp(v, up);
      double sum = 0;
      double count = 0;
      for (int i = 0; i < up.n; i++) {
	ma::Entity* currentEdge = up.e[i];
	ma::Entity* vs[2];
	mesh->getDownward(currentEdge, 0, vs);
	ma::Vector diff = ma::getPosition(mesh, vs[0]) - ma::getPosition(mesh, vs[1]);
	sum += diff.getLength();
	count++;
      }
      double currentSize = sum/count;

      // compute the size multiplier
      double y = 1./2.;
      int tag  = mesh->getModelTag(mesh->toModel(v));
      int type = mesh->getModelType(mesh->toModel(v));
      if (type == 2 && (tag == 66 || tag == 86))
      	y = 3;
      if (type == 1 && tag == 178)
      	y = 3;
      // new size is multiplier * currentSize
      return y*currentSize;
    }
  private:
    ma::Mesh* mesh;
    double factor;
    double average;
    ma::Vector upper;
    ma::Vector lower;
};


class GeomRobin : public ma::IsotropicFunction
{
  public:
    GeomRobin(ma::Mesh* m, double inFactor)
    {
      mesh = m;
      factor = inFactor;
      average = ma::getAverageEdgeLength(m);
      ma::getBoundingBox(m, lower, upper);
      printf("average is %f\n", average);
    }
    virtual double getValue(ma::Entity* v)
    {
      ma::Vector p = ma::getPosition(mesh, v);
      double x = p[0];
      (void)x;

      // get current size
      apf::Up up;
      mesh->getUp(v, up);
      double sum = 0;
      double count = 0;
      for (int i = 0; i < up.n; i++) {
	ma::Entity* currentEdge = up.e[i];
	ma::Entity* vs[2];
	mesh->getDownward(currentEdge, 0, vs);
	ma::Vector diff = ma::getPosition(mesh, vs[0]) - ma::getPosition(mesh, vs[1]);
	sum += diff.getLength();
	count++;
      }
      double currentSize = sum/count;

      // compute the size multiplier
      double y = 1./2.;
      int tag  = mesh->getModelTag(mesh->toModel(v));
      int type = mesh->getModelType(mesh->toModel(v));
      if (type == 2 && (tag == 13 || tag == 26))
      	y = 3.;
      /* if (type == 1 && tag == 178) */
      /* 	y = 3; */
      // new size is multiplier * currentSize
      return y*currentSize;
    }
  private:
    ma::Mesh* mesh;
    double factor;
    double average;
    ma::Vector upper;
    ma::Vector lower;
};

