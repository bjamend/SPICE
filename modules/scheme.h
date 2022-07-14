// -----------------------------------------------------------------------------
// Evolve the simulation in time (RK scheme, 3rd order).
// -----------------------------------------------------------------------------


#include "boundaries.h"
#include "riemann.h"


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
void rk3(double *u0, double t, double *x, double *y, double dx, double dy, double *events, int counter) {

  while (t < t_final) {

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
                                  counter) * dt;
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
                                                counter) * dt;
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
                                                counter) * dt;
      }
    }
    free(u2);
    t = (1.0 / 3.0) * t0 + (2.0 / 3.0) * (t + dt);
    memcpy(u0, u3, num_zones * num_zones * sizeof(double));
    free(u3);

  printf("t=%f\n", t);
  counter += 1;

  }
}
