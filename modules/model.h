// -----------------------------------------------------------------------------
// Prescriptions for diffusion coefficients and flow velocities, as well as
// an initial mass concentration.
// -----------------------------------------------------------------------------


// spatially-inhomogeneous large-scale flow velocity
double flow_velocity(double x, double y, char axis) {
  double vx = 0.0;
  double vy = 0.0;
  if (axis == 'x') {
    return vx;
  } if (axis == 'y') {
    return vy;
  }
  return 0.0;
}


// spatially-inhomogeneous diffusion coefficient
double diffusion_coefficient(double x, double y) {
  return 0.2;
}


// establish initial u (should be 0 by default)
double initial_condition(double x, double y) {
  // return sin(10.0 * 3.14159 * x) * sin(10.0 * 3.14159 * y);
  return 0.0;
}
