// -----------------------------------------------------------------------------
// Construct a grid of vertices.
// -----------------------------------------------------------------------------


#include "../input.h"


// populates arrays with vertex coordinates
void construct_grid(double *x, double *y, double x_l, double y_l,
                    double dx, double dy){
  for (int i = 0; i < num_zones; ++i) {
    x[i] = x_l + (i + 0.5) * dx;
    y[i] = y_l + (i + 0.5) * dy;
  }
}
