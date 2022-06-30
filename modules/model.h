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
  return 0.05;
}


// establish initial u
double initial_condition(double x, double y) {
  if (sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)) < 0.025) {
    return 1.0;
  }
  return 0.0;
}
