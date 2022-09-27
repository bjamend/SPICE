#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


// -----------------------------------------------------------------------------
// Various types of boundary conditions.
// -----------------------------------------------------------------------------


void boundary_condition(double *u0, double *u1, double *u2, double *u3,
                        int num_zones, char type) {
  // f = Fixed
  // o = Outflow

  if (type == 'f') {
    for (int i = 0; i < num_zones; ++i) {
      u1[0 * num_zones + i] = u0[0 * num_zones + i];
      u1[1 * num_zones + i] = u0[1 * num_zones + i];
      u1[(num_zones - 1) * num_zones + i] = u0[(num_zones - 1) * num_zones + i];
      u1[(num_zones - 2) * num_zones + i] = u0[(num_zones - 2) * num_zones + i];
      u1[i * num_zones + 0] = u0[i * num_zones + 0];
      u1[i * num_zones + 1] = u0[i * num_zones + 1];
      u1[i * num_zones + (num_zones - 1)] = u0[i * num_zones + (num_zones - 1)];
      u1[i * num_zones + (num_zones - 2)] = u0[i * num_zones + (num_zones - 2)];

      u2[0 * num_zones + i] = u0[0 * num_zones + i];
      u2[1 * num_zones + i] = u0[1 * num_zones + i];
      u2[(num_zones - 1) * num_zones + i] = u0[(num_zones - 1) * num_zones + i];
      u2[(num_zones - 2) * num_zones + i] = u0[(num_zones - 2) * num_zones + i];
      u2[i * num_zones + 0] = u0[i * num_zones + 0];
      u2[i * num_zones + 1] = u0[i * num_zones + 1];
      u2[i * num_zones + (num_zones - 1)] = u0[i * num_zones + (num_zones - 1)];
      u2[i * num_zones + (num_zones - 2)] = u0[i * num_zones + (num_zones - 2)];

      u3[0 * num_zones + i] = u0[0 * num_zones + i];
      u3[1 * num_zones + i] = u0[1 * num_zones + i];
      u3[(num_zones - 1) * num_zones + i] = u0[(num_zones - 1) * num_zones + i];
      u3[(num_zones - 2) * num_zones + i] = u0[(num_zones - 2) * num_zones + i];
      u3[i * num_zones + 0] = u0[i * num_zones + 0];
      u3[i * num_zones + 1] = u0[i * num_zones + 1];
      u3[i * num_zones + (num_zones - 1)] = u0[i * num_zones + (num_zones - 1)];
      u3[i * num_zones + (num_zones - 2)] = u0[i * num_zones + (num_zones - 2)];
    }
  }
  if (type == 'o') {
    for (int i = 0; i < num_zones; ++i) {
      u1[0 * num_zones + i] = u0[2 * num_zones + i];
      u1[1 * num_zones + i] = u0[2 * num_zones + i];
      u1[(num_zones - 1) * num_zones + i] = u0[(num_zones - 3) * num_zones + i];
      u1[(num_zones - 2) * num_zones + i] = u0[(num_zones - 3) * num_zones + i];
      u1[i * num_zones + 0] = u0[i * num_zones + 2];
      u1[i * num_zones + 1] = u0[i * num_zones + 2];
      u1[i * num_zones + (num_zones - 1)] = u0[i * num_zones + (num_zones - 3)];
      u1[i * num_zones + (num_zones - 2)] = u0[i * num_zones + (num_zones - 3)];

      u2[0 * num_zones + i] = u0[2 * num_zones + i];
      u2[1 * num_zones + i] = u0[2 * num_zones + i];
      u2[(num_zones - 1) * num_zones + i] = u0[(num_zones - 3) * num_zones + i];
      u2[(num_zones - 2) * num_zones + i] = u0[(num_zones - 3) * num_zones + i];
      u2[i * num_zones + 0] = u0[i * num_zones + 2];
      u2[i * num_zones + 1] = u0[i * num_zones + 2];
      u2[i * num_zones + (num_zones - 1)] = u0[i * num_zones + (num_zones - 3)];
      u2[i * num_zones + (num_zones - 2)] = u0[i * num_zones + (num_zones - 3)];

      u3[0 * num_zones + i] = u0[2 * num_zones + i];
      u3[1 * num_zones + i] = u0[2 * num_zones + i];
      u3[(num_zones - 1) * num_zones + i] = u0[(num_zones - 3) * num_zones + i];
      u3[(num_zones - 2) * num_zones + i] = u0[(num_zones - 3) * num_zones + i];
      u3[i * num_zones + 0] = u0[i * num_zones + 2];
      u3[i * num_zones + 1] = u0[i * num_zones + 2];
      u3[i * num_zones + (num_zones - 1)] = u0[i * num_zones + (num_zones - 3)];
      u3[i * num_zones + (num_zones - 2)] = u0[i * num_zones + (num_zones - 3)];
    }
  }
}


// -----------------------------------------------------------------------------
// Additional math functions not included in standard libraries.
// -----------------------------------------------------------------------------


// returns the sign of a value
double sign(double x) {
  if (x < 0.0) {
    return -1.0;
  } else if (x > 0.0) {
    return 1.0;
  } else {
    return 0.0;
  }
}


// returns the minimum of three values
double min3(double a, double b, double c) {
  if (a < b) {
    if (a < c) {
      return a;
    } else if (a > c) {
      return c;
    } else {
      return a;
    }
  } else if (b < a) {
    if (b < c) {
      return b;
    } else if (b > c) {
      return c;
    } else {
      return c;
    }
  } else if (a == b) {
    if (c < b) {
      return c;
    } else {
      return b;
    };
  } else {
    return 0.0;
  }
}


// slope limiting function
double minmod(double a, double b, double c) {
  return 0.25 * fabs(sign(a) + sign(b)) * (sign(a) + sign(c)) *
         min3(fabs(a), fabs(b), fabs(c));
}


// computes a slope-limited gradient
double plm_gradient(double a, double b, double c) {
  double plm_theta = 1.5;
  return minmod(plm_theta * (b - a), 0.5 * (c - a), plm_theta * (c - b));
}


// generate random number
double random_number(double min, double max) {
  double range = (max - min);
  double div   = RAND_MAX / range;
  return min + (rand() / div);
}


// -----------------------------------------------------------------------------
// Construct a grid of vertices.
// -----------------------------------------------------------------------------


#include "input.h"


// populates arrays with vertex coordinates
void construct_grid(double *x, double *y, double x_l, double y_l,
                    double dx, double dy){
  for (int i = 0; i < num_zones; ++i) {
    x[i] = x_l + (i + 0.5) * dx;
    y[i] = y_l + (i + 0.5) * dy;
  }
}


// -----------------------------------------------------------------------------
// Save data to text file.
// -----------------------------------------------------------------------------


// exports coordinates, initial u, and final u to a text file
void export_data(double *x, double *y, double *u, double *u0, int num_zones,
                 int counter) {
  char filepath[256];
  snprintf(filepath, sizeof(filepath), "data/data%d.txt", counter);
  FILE *f = fopen(filepath, "w");
  if (f == NULL) {
    printf("Error opening file!\n");
    exit(1);
  }
  for (int i = 0; i < num_zones; ++i) {
    for (int j = 0; j < num_zones; ++j) {
      fprintf(f, "%f, %f, %f, %f\n", x[i], y[j], u[i*num_zones+j],
              u0[i*num_zones+j]);
    }
  }
  fclose(f);
}

// import spatial and temporal coordinates for r-process events
void generate_array(double *sites) {
  FILE *file;
  file = fopen("data/eventsites.csv", "r");
  double num;
  int i = 0;
  while(fscanf(file, "%lf", &num) == 1) {
    sites[i] = num;
    i++;
  }
  fclose(file);
}


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


// -----------------------------------------------------------------------------
// Compute fluxes to solve the Riemann problem (2nd-order).
// -----------------------------------------------------------------------------


// flux from advective term in diffusive-advective equation
double advective_flux(double ul, double ur, double x, double y, char axis) {
  double vx = flow_velocity(x, y, 'x');
  double vy = flow_velocity(x, y, 'y');
  if (axis == 'x') {
    if (vx < 0.0) {
      return vx * ur;
    } if (vx > 0.0) {
      return vx * ul;
    }
  } if (axis == 'y') {
    if (vy < 0.0) {
      return vy * ur;
    } if (vy > 0.0) {
      return vy * ul;
    }
  }
  return 0.0;
}


// source terms
double source_terms(double x, double y, double t, double *events, int counter0,
                    int counter) {
  double result = 0.0;
  for (int i = counter0; i < counter; ++i) {
    double x_0 = events[3*i+1];
    double y_0 = events[3*i+2];
    if (((x - x_0) * (x - x_0) + (y - y_0) * (y - y_0)) < (25.0 * r * r)) {
      result += 10000.0 * exp(-((x - x_0) * (x - x_0) + (y - y_0) * (y - y_0))
                              / r / r);
    }
  }
  return result;
}


// flux from diffusive term in diffusive-advective equation
double diffusive_flux(double ul, double ur, double x, double y, double dx,
                      double dy, char axis) {
  if (axis == 'x') {
    return -diffusion_coefficient(x, y) * (ur - ul) / dx;
  } if (axis == 'y') {
    return -diffusion_coefficient(x, y) * (ur - ul) / dy;
  }
  return 0.0;
}


// time derivative of u
double du_dt(double u_im2j, double u_im1j, double u_ij, double u_ip1j,
             double u_ip2j, double u_ijm2, double u_ijm1, double u_ijp1,
             double u_ijp2, double x, double y, double dx, double dy,
             double t, double *events, int counter0, int counter) {

  double u_limhj = u_im1j + 0.5 * plm_gradient(u_im2j, u_im1j, u_ij);
  double u_rimhj = u_ij   - 0.5 * plm_gradient(u_im1j, u_ij, u_ip1j);
  double u_liphj = u_ij   + 0.5 * plm_gradient(u_im1j, u_ij, u_ip1j);
  double u_riphj = u_ip1j - 0.5 * plm_gradient(u_ij, u_ip1j, u_ip2j);
  double u_lijmh = u_ijm1 + 0.5 * plm_gradient(u_ijm2, u_ijm1, u_ij);
  double u_rijmh = u_ij   - 0.5 * plm_gradient(u_ijm1, u_ij, u_ijp1);
  double u_lijph = u_ij   + 0.5 * plm_gradient(u_ijm1, u_ij, u_ijp1);
  double u_rijph = u_ijp1 - 0.5 * plm_gradient(u_ij, u_ijp1, u_ijp2);


  double f_iphj = advective_flux(u_liphj, u_riphj, x + 0.5 * dx, y, 'x') +
                  diffusive_flux(u_ij, u_ip1j, x + 0.5 * dx, y, dx, dy, 'x');
  double f_imhj = advective_flux(u_limhj, u_rimhj, x - 0.5 * dx, y, 'x') +
                  diffusive_flux(u_im1j, u_ij, x - 0.5 * dx, y, dx, dy, 'x');
  double g_ijph = advective_flux(u_lijph, u_rijph, x, y + 0.5 * dy, 'y') +
                  diffusive_flux(u_ij, u_ijp1, x, y + 0.5 * dy, dx, dy, 'y');
  double g_ijmh = advective_flux(u_lijmh, u_rijmh, x, y - 0.5 * dy, 'y') +
                  diffusive_flux(u_ijm1, u_ij, x, y - 0.5 * dy, dx, dy, 'y');
  double source = source_terms(x, y, t, events, counter0, counter);
  return -(f_iphj - f_imhj) / dx - (g_ijph - g_ijmh) / dy + source;
}


// -----------------------------------------------------------------------------
// Evolve the simulation in time (RK scheme, 3rd order).
// -----------------------------------------------------------------------------


// minimum timestep based on flow velocity
double advective_timestep(double *x, double *y, double dx, double dy) {
  double max_velocity = 0.0;
  for (int i = 0; i < num_zones; ++i) {
    for (int j = 0; j < num_zones; ++j) {
      double net_flow_velocity = sqrt(pow(flow_velocity(x[i], y[j], 'x'), 2.0) +
                                      pow(flow_velocity(x[i], y[j], 'y'), 2.0));
      if (net_flow_velocity > max_velocity) {
        max_velocity = net_flow_velocity;
      }
    }
  }
  return sqrt(dx * dy) / (max_velocity + 1.0E-12);
}


//minimum timestep based on diffusion coefficient
double diffusive_timestep(double *x, double *y, double dx, double dy) {
  double max_diffusion = 0.0;
  for (int i = 0; i < num_zones; ++i) {
    for (int j = 0; j < num_zones; ++j){
      if (diffusion_coefficient(x[i], y[j]) > max_diffusion) {
        max_diffusion = diffusion_coefficient(x[i], y[j]);
      }
    }
  }
  return 0.5 * (dx * dy) / (max_diffusion + 1.0E-12);
}


// overall minimum timestep
double timestep(double *x, double *y, double dx, double dy) {
  double a_t = advective_timestep(x, y, dx, dy);
  double d_t = diffusive_timestep(x, y, dx, dy);
  if (a_t < d_t) {
    return a_t;
  }
  return d_t;
}


// main rk3 algorithm
void rk3(double *u0, double t, double *x, double *y, double dx, double dy,
         double *events) {

  int counter = 0;
  int time_counter = 0;

  while (t < t_final) {

    double counter0 = counter;

    while (events[3*counter] < t) {
      counter += 1;
    }

    // printf("%f\n", counter-counter0);

    double dt  = cfl * timestep(x, y, dx, dy);
    double *u1 = malloc(num_zones * num_zones * sizeof(double));
    double *u2 = malloc(num_zones * num_zones * sizeof(double));
    double *u3 = malloc(num_zones * num_zones * sizeof(double));
    double t0  = t;

    boundary_condition(u0, u1, u2, u3, num_zones, 'o');

    for (int i = 2; i < (num_zones-2); ++i) {
      for (int j = 2; j < (num_zones-2); ++j) {
        u1[i*num_zones+j] = u0[i*num_zones+j] +
                            du_dt(u0[(i-2)*num_zones+(j+0)],
                                  u0[(i-1)*num_zones+(j+0)],
                                  u0[(i+0)*num_zones+(j+0)],
                                  u0[(i+1)*num_zones+(j+0)],
                                  u0[(i+2)*num_zones+(j+0)],
                                  u0[(i+0)*num_zones+(j-2)],
                                  u0[(i+0)*num_zones+(j-1)],
                                  u0[(i+0)*num_zones+(j+1)],
                                  u0[(i+0)*num_zones+(j+2)],
                                  x[i], y[j], dx, dy, t, events,
                                  counter0, counter) * dt;
      }
    }
    t += dt;

    for (int i = 2; i < (num_zones-2); ++i) {
      for (int j = 2; j < (num_zones-2); ++j) {
        u2[i*num_zones+j] = (3.0 / 4.0) * u0[i*num_zones+j] +
                            (1.0 / 4.0) * u1[i*num_zones+j] +
                            (1.0 / 4.0) * du_dt(u1[(i-2)*num_zones+(j+0)],
                                                u1[(i-1)*num_zones+(j+0)],
                                                u1[(i+0)*num_zones+(j+0)],
                                                u1[(i+1)*num_zones+(j+0)],
                                                u1[(i+2)*num_zones+(j+0)],
                                                u1[(i+0)*num_zones+(j-2)],
                                                u1[(i+0)*num_zones+(j-1)],
                                                u1[(i+0)*num_zones+(j+1)],
                                                u1[(i+0)*num_zones+(j+2)],
                                                x[i], y[j], dx, dy, t, events,
                                                counter0, counter) * dt;
      }
    }
    free(u1);
    t = (3.0 / 4.0) * t0 + (1.0 / 4.0) * (t + dt);

    for (int i = 2; i < (num_zones-2); ++i) {
      for (int j = 2; j < (num_zones-2); ++j) {
        u3[i*num_zones+j] = (1.0 / 3.0) * u0[i*num_zones+j] +
                            (2.0 / 3.0) * u2[i*num_zones+j] +
                            (2.0 / 3.0) * du_dt(u2[(i-2)*num_zones+(j+0)],
                                                u2[(i-1)*num_zones+(j+0)],
                                                u2[(i+0)*num_zones+(j+0)],
                                                u2[(i+1)*num_zones+(j+0)],
                                                u2[(i+2)*num_zones+(j+0)],
                                                u2[(i+0)*num_zones+(j-2)],
                                                u2[(i+0)*num_zones+(j-1)],
                                                u2[(i+0)*num_zones+(j+1)],
                                                u2[(i+0)*num_zones+(j+2)],
                                                x[i], y[j], dx, dy, t, events,
                                                counter0, counter) * dt;
      }
    }
    free(u2);
    t = (1.0 / 3.0) * t0 + (2.0 / 3.0) * (t + dt);
    memcpy(u0, u3, num_zones * num_zones * sizeof(double));
    free(u3);

    if ((time_counter % 100) == 0) {
      export_data(x, y, u0, u0, num_zones, time_counter);
    }

    time_counter += 1;

  printf("t=%f, %d\n", t, time_counter);

  }
}


int main() {

  // initialize spatial grid, time, and concentration
  double *events = malloc(3 * 10000 * sizeof(double));
  generate_array(events);
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
  rk3(u0, t, x, y, dx, dy, events);

  free(u0);
  free(ui);

  return 0;
}
