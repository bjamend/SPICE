// -----------------------------------------------------------------------------
// Construct a grid of vertices.
// -----------------------------------------------------------------------------


#include "../input.h"


const double x_l     = 0.0;
const double x_r     = 1.0;
const double dx      = (x_r - x_l) / num_zones;
const double y_l     = 0.0;
const double y_r     = 1.0;
const double dy      = (y_r - y_l) / num_zones;


// populates arrays with vertex coordinates
void construct_grid(double *x, double *y){
  for (int i = 0; i < num_zones; ++i) {
    x[i] = x_l + (i + 0.5) * dx;
    y[i] = y_l + (i + 0.5) * dy;
  }
}
