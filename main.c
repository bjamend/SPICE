#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "modules/io.h"
#include "modules/grid.h"
#include "modules/scheme.h"


int main() {

  // initialize spatial grid, time, and concentration
  double x[num_zones];
  double y[num_zones];
  const double x_l     = 0.0;
  const double x_r     = 1.0;
  const double dx      = (x_r - x_l) / num_zones;
  const double y_l     = 0.0;
  const double y_r     = 1.0;
  const double dy      = (y_r - y_l) / num_zones;
  construct_grid(x, y, x_l, y_l, dx, dy);
  double t   = 0.0;
  double *u0 = malloc(num_zones * num_zones * sizeof(double));
  double *ui = malloc(num_zones * num_zones * sizeof(double));
  for (int i = 0; i < num_zones; ++i) {
    for (int j = 0; j < num_zones; ++j) {
      u0[i*num_zones+j] = initial_condition(x[i], y[j]);
      ui[i*num_zones+j] = initial_condition(x[i], y[j]);
    }
  }

  // evolve the simulation in time
  rk3(u0, t, x, y, dx, dy);

  // store initial and final data in a text file
  export_data(x, y, u0, ui, num_zones);

  free(u0);
  free(ui);
  return 0;
}
