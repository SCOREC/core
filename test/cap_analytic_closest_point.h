void doubleConeClosestPointAnalytic(double const from[3], double to[3], double to_norm[3]) {
  double x0 = from[0];
  // x axis axisymmetry, y0 here is distance from x axis
  double y0 = std::sqrt(from[1]*from[1] + from[2]*from[2]);

  // For Mathematica CForm output
  #define Power(base, exp) std::pow(base,exp)
  #define Sqrt(arg) std::sqrt(arg)

  // from NX, intersection point Point( -3.27529516[m], 1.84356815[m], 0.0[m] )
  double intersection_y = 1.84356815;
  // Parabola
  double pA =  -0.0516559000000000;
  double pZ = -3.0997300000000000;
  double p_y = (-1 - 2*pA*pZ + 2*pA*x0)/(Power(6,0.3333333333333333)*Power(9*Power(pA,4)*y0 + Sqrt(3)*Sqrt(Power(pA,6)*(2*Power(1 + 2*pA*(pZ - x0),3) + 27*Power(pA,2)*Power(y0,2))),0.3333333333333333)) + Power(9*Power(pA,4)*y0 + Sqrt(3)*Sqrt(Power(pA,6)*(2*Power(1 + 2*pA*(pZ - x0),3) + 27*Power(pA,2)*Power(y0,2))),0.3333333333333333)/ (Power(6,0.6666666666666666)*Power(pA,2));
  p_y = std::max(intersection_y, p_y);
  double p_x = pA*p_y*p_y + pZ;
  double p_d2 = std::pow(p_x-x0,2) + std::pow(p_y-y0,2);

  // Cone (l for line)
  double A = 1000;
  double B = 1786.99;
  double C = -19.1427;
  //double l_d2 = Power(C + A*x0 + B*y0,2)/(Power(A,2) + Power(B,2));
  //double l_x = (-(A*C) + B*(B*x0 - A*y0))/(Power(A,2) + Power(B,2));
  double l_y = (-(B*C) + A*(-(B*x0) + A*y0))/(Power(A,2) + Power(B,2));
  l_y = std::min(intersection_y, l_y);
  l_y = std::max(0.0d, l_y);
  double l_x = -(B*l_y + C)/A;
  double l_d2 = std::pow(l_x-x0,2) + std::pow(l_y-y0,2);

  /* old idea for rejecting invalid cone/parabola regions
  p_x_l_y = pA*l_y*l_y + pZ; // x coordinate of closest point on cone projected in the x direction to the parabola
  l_x_p_y = -(B*p_y + C)/A; // x coordinate of the closest point on parabola projected in the x direction to the cone
  if (p_y > l_y) {
    // outside cone-parabola intersection
    // if the closest point on the cone is behind the
    // closest point on the parabola in this region
    // make parabola win
    if (l_x_p_y < p_x) p_d2 = 0; 
  } else {
    // inside cone-parabola intersection
    // if closest point on the parabola is behind the
    // cloest point on the cone in this region
    // make cone win
    if (p_x_l_y < l_x) l_d2 = 0;
  }
  */

  // Convert to x y z
  // Normals
  double cls_y;
  double dx;
  double dy;
  if (p_d2 < l_d2) {
    to[0] = p_x;
    cls_y = p_y; 
    dx = -1.0;
    dy = 2*pA*p_y;
  } else {
    to[0] = l_x;
    cls_y = l_y;
    dx = A;
    dy = B;
  }

  double zero_tol = 1e-6;
  double ratio = cls_y/std::max(y0,zero_tol);
  to[1] = from[1]*ratio;
  to[2] = from[2]*ratio;

  // Build normal vector in x y z
  to_norm[0] = dx;
  double norm_ratio = dy/std::max(y0,zero_tol);
  to_norm[1] = from[1]*norm_ratio;
  to_norm[2] = from[2]*norm_ratio;
  double norm_norm = std::sqrt(to_norm[0]*to_norm[0] + to_norm[1]*to_norm[1] + to_norm[2]*to_norm[2]);
  to_norm[0] /= norm_norm;
  to_norm[1] /= norm_norm;
  to_norm[2] /= norm_norm;
}