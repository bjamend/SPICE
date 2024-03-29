#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "input.h"

#define sign(x)       ((x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0))
#define min2(a, b)    (a < b ? a : b)
#define min3(a, b, c) (min2(a, min2(b, c)))

#define MIN_TIMESTEP_THRESHOLD 1E-12
#define PI 3.14159


struct Source {
    double x, y;
};


void boundary_condition(double *u0, double *u1, double *u2, double *u3,
                        int num_zones)
{
    for (int i = 0; i < num_zones; ++i) {
      for (int l = 0; l < 2; ++l) {
        u1[i * num_zones + l] = u0[i * num_zones + num_zones - 4 + l];
        u1[i * num_zones + num_zones - l - 1] = u0[i * num_zones + 3 + l];

        u2[i * num_zones + l] = u0[i * num_zones + num_zones - 4 + l];
        u2[i * num_zones + num_zones - l - 1] = u0[i * num_zones + 3 + l];

        u3[i * num_zones + l] = u0[i * num_zones + num_zones - 4 + l];
        u3[i * num_zones + num_zones - l - 1] = u0[i * num_zones + 3 + l];
      }
    }

    for (int j = 0; j < num_zones; ++j) {
      for (int l = 0; l < 2; ++l) {
        u1[l * num_zones + j] = u0[(num_zones - 4 + l) * num_zones + j];
        u1[(num_zones - l - 1) * num_zones + j] = u0[(3 - l) * num_zones + j];

        u2[l * num_zones + j] = u0[(num_zones - 4 + l) * num_zones + j];
        u2[(num_zones - l - 1) * num_zones + j] = u0[(3 - l) * num_zones + j];

        u3[l * num_zones + j] = u0[(num_zones - 4 + l) * num_zones + j];
        u3[(num_zones - l - 1) * num_zones + j] = u0[(3 - l) * num_zones + j];
      }
    }
}

// slope limiting function
double minmod(double a, double b, double c)
{
    return 0.25 * fabs(sign(a) + sign(b)) * (sign(a) + sign(c)) *
           min3(fabs(a), fabs(b), fabs(c));
}

// computes a slope-limited gradient
double plm_gradient(double a, double b, double c)
{
    double plm_theta = 1.5;
    return minmod(plm_theta * (b - a), 0.5 * (c - a), plm_theta * (c - b));
}

// computes the factorial of an integer n
int factorial(int n)
{
    if (n == 0) {
        return 1;
    }
    return n * factorial(n - 1);
}

// generates a random double between 0 and 1
double random_double() { return (double)rand() / (double)RAND_MAX; }

// generates samples from a poisson distribution
int poisson_sample(double lambda, int max_count)
{
    double u = random_double();
    double p = 0;
    int count = 0;

    while (count < max_count) {
        p += (pow(lambda, count) * exp(-lambda) / factorial(count));
        if (u <= p) {
            return count;
        }
        count += 1;
    }
    return 0;
}

// populates arrays with vertex coordinates
void construct_grid(double *x, double *y, double x_l, double y_l, double dx,
                    double dy)
{
    for (int i = 0; i < num_zones; ++i) {
        x[i] = x_l + (i + 0.5) * dx;
    }
    for (int j = 0; j < num_zones; ++j) {
        y[j] = y_l + (j + 0.5) * dy;
    }
}

// exports coordinates, initial u, and final u to a text file
void export_data(double t, double *x, double *y, double *u, int num_zones,
                 int counter)
{
    char filepath[256];
    snprintf(filepath, sizeof(filepath), "data/checkpoint%d.txt", counter / checkpoint_interval);
    FILE *f = fopen(filepath, "w");
    if (f == NULL) {
        printf("Error opening file (did you create a 'data' directory?)\n");
        exit(1);
    }
    for (int i = 0; i < num_zones; ++i) {
        for (int j = 0; j < num_zones; ++j) {
            fprintf(f, "%.10e, %f, %f, %.10e\n", t, x[i], y[j],
                    u[i * num_zones + j]);
        }
    }
    fclose(f);
}

// spatially-inhomogeneous large-scale flow velocity
double flow_velocity(double x, double y, char axis)
{
    double vx = 0.0;
    double vy = 0.0;

    if (axis == 'x') {
        return vx;
    }
    if (axis == 'y') {
        return vy;
    }
    return 0.0;
}

// spatially-inhomogeneous diffusion coefficient
double diffusion_coefficient(double x, double y) { return 0.2; }

// establish initial u (should be 0 by default)
double initial_condition(double x, double y)
{
  return 0.0;
}

// flux from advective term in diffusive-advective equation
double advective_flux(double ul, double ur, double x, double y, char axis)
{
    double vx = flow_velocity(x, y, 'x');
    double vy = flow_velocity(x, y, 'y');
    if (axis == 'x') {
        if (vx < 0.0) {
            return vx * ur;
        }
        if (vx > 0.0) {
            return vx * ul;
        }
    }
    if (axis == 'y') {
        if (vy < 0.0) {
            return vy * ur;
        }
        if (vy > 0.0) {
            return vy * ul;
        }
    }
    return 0.0;
}

// flux from diffusive term in diffusive-advective equation
double diffusive_flux(double ul, double ur, double x, double y, double dx,
                      double dy, char axis)
{
    if (axis == 'x') {
        return -diffusion_coefficient(x, y) * (ur - ul) / dx;
    }
    if (axis == 'y') {
        return -diffusion_coefficient(x, y) * (ur - ul) / dy;
    }
    return 0.0;
}

double source_kernel(double x, double y, double r, double x_0, double y_0, double dt) {
  double radius = pow((x - x_0) * (x - x_0) + (y - y_0) * (y - y_0), 0.5);
  if (radius < r) {
    return source_yield / PI / r / r / dt;
  }
  return 0.0;
}

void generate_source_coordinates(struct Source *sources,
                                 int *num_sources_this_step, double dt)
{
    double lambda = dt * source_production_rate;
    int source_count = poisson_sample(lambda, max_sources_per_timestep);

    for (int i = 0; i < source_count; ++i) {
      sources[i].x = 50.0 * random_double();
      sources[i].y = 50.0 * random_double();
    }

    *num_sources_this_step = source_count;
}

// time derivative of u
double du_dt(double u_im2j, double u_im1j, double u_ij, double u_ip1j,
             double u_ip2j, double u_ijm2, double u_ijm1, double u_ijp1,
             double u_ijp2, double x, double y, double dx, double dy, double t,
             double dt, struct Source *sources, int num_sources)
{
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

    double source_terms = 0.0;
    for (int i = 0; i < num_sources; ++i) {
      double x_0 = sources[i].x;
      double y_0 = sources[i].y;
      source_terms += source_kernel(x, y, r, x_0, y_0, dt);
    }

    return -(f_iphj - f_imhj) / dx - (g_ijph - g_ijmh) / dy + source_terms;
}

// minimum timestep based on flow velocity
double advective_timestep(double *x, double *y, double dx, double dy)
{
    double max_velocity = 0.0;
    for (int i = 0; i < num_zones; ++i) {
        for (int j = 0; j < num_zones; ++j) {
            double net_flow_velocity =
                sqrt(pow(flow_velocity(x[i], y[j], 'x'), 2.0) +
                     pow(flow_velocity(x[i], y[j], 'y'), 2.0));
            if (net_flow_velocity > max_velocity) {
                max_velocity = net_flow_velocity;
            }
        }
    }
    return sqrt(dx * dy) / (max_velocity + MIN_TIMESTEP_THRESHOLD);
}

// minimum timestep based on diffusion coefficient
double diffusive_timestep(double *x, double *y, double dx, double dy)
{
    double max_diffusion = 0.0;
    for (int i = 0; i < num_zones; ++i) {
        for (int j = 0; j < num_zones; ++j) {
            if (diffusion_coefficient(x[i], y[j]) > max_diffusion) {
                max_diffusion = diffusion_coefficient(x[i], y[j]);
            }
        }
    }
    return 0.05 * (dx * dy) / (max_diffusion + MIN_TIMESTEP_THRESHOLD);
}

// overall minimum timestep
double timestep(double *x, double *y, double dx, double dy)
{
    double a_t = advective_timestep(x, y, dx, dy);
    double d_t = diffusive_timestep(x, y, dx, dy);
    if (a_t < d_t) {
        return a_t;
    }
    return d_t;
}

void abort_if_nan(double *u) {
    for (int i = 0; i < (num_zones * num_zones); ++i) {
      if (isnan(u[i])) {
        printf("Code crashed, force exit.\n");
        exit(1);
      }
    }
}

// main rk3 algorithm
void rk3(double *u0, double t, double *x, double *y, double dx, double dy)
{
    int time_counter = 0;

    struct Source *sources = (struct Source *)malloc(max_sources_per_timestep *
                                                     sizeof(struct Source));
    int num_sources_this_step = 0;

    double *u1 = malloc(num_zones * num_zones * sizeof(double));
    double *u2 = malloc(num_zones * num_zones * sizeof(double));
    double *u3 = malloc(num_zones * num_zones * sizeof(double));

    while (t < t_final) {

        clock_t start = clock();

        double dt = cfl * timestep(x, y, dx, dy);
        double t0 = t;

        boundary_condition(u0, u1, u2, u3, num_zones);
        generate_source_coordinates(sources, &num_sources_this_step, dt);

        for (int i = 2; i < (num_zones - 2); ++i) {
            for (int j = 2; j < (num_zones - 2); ++j) {
                u1[i * num_zones + j] =
                    u0[i * num_zones + j] +
                    du_dt(u0[(i - 2) * num_zones + (j + 0)],
                          u0[(i - 1) * num_zones + (j + 0)],
                          u0[(i + 0) * num_zones + (j + 0)],
                          u0[(i + 1) * num_zones + (j + 0)],
                          u0[(i + 2) * num_zones + (j + 0)],
                          u0[(i + 0) * num_zones + (j - 2)],
                          u0[(i + 0) * num_zones + (j - 1)],
                          u0[(i + 0) * num_zones + (j + 1)],
                          u0[(i + 0) * num_zones + (j + 2)], x[i], y[j], dx, dy,
                          t, dt, sources, num_sources_this_step) * dt;
            }
        }
        t += dt;

        for (int i = 2; i < (num_zones - 2); ++i) {
            for (int j = 2; j < (num_zones - 2); ++j) {
                u2[i * num_zones + j] =
                    (3.0 / 4.0) * u0[i * num_zones + j] +
                    (1.0 / 4.0) * u1[i * num_zones + j] +
                    (1.0 / 4.0) *
                        du_dt(u1[(i - 2) * num_zones + (j + 0)],
                              u1[(i - 1) * num_zones + (j + 0)],
                              u1[(i + 0) * num_zones + (j + 0)],
                              u1[(i + 1) * num_zones + (j + 0)],
                              u1[(i + 2) * num_zones + (j + 0)],
                              u1[(i + 0) * num_zones + (j - 2)],
                              u1[(i + 0) * num_zones + (j - 1)],
                              u1[(i + 0) * num_zones + (j + 1)],
                              u1[(i + 0) * num_zones + (j + 2)], x[i], y[j], dx,
                              dy, t, dt, sources, num_sources_this_step) * dt;
            }
        }

        t = (3.0 / 4.0) * t0 + (1.0 / 4.0) * (t + dt);

        for (int i = 2; i < (num_zones - 2); ++i) {
            for (int j = 2; j < (num_zones - 2); ++j) {
                u3[i * num_zones + j] =
                    (1.0 / 3.0) * u0[i * num_zones + j] +
                    (2.0 / 3.0) * u2[i * num_zones + j] +
                    (2.0 / 3.0) *
                        du_dt(u2[(i - 2) * num_zones + (j + 0)],
                              u2[(i - 1) * num_zones + (j + 0)],
                              u2[(i + 0) * num_zones + (j + 0)],
                              u2[(i + 1) * num_zones + (j + 0)],
                              u2[(i + 2) * num_zones + (j + 0)],
                              u2[(i + 0) * num_zones + (j - 2)],
                              u2[(i + 0) * num_zones + (j - 1)],
                              u2[(i + 0) * num_zones + (j + 1)],
                              u2[(i + 0) * num_zones + (j + 2)], x[i], y[j], dx,
                              dy, t, dt, sources, num_sources_this_step) * dt;
            }
        }
        t = (1.0 / 3.0) * t0 + (2.0 / 3.0) * (t + dt);
        memcpy(u0, u3, num_zones * num_zones * sizeof(double));

        abort_if_nan(u0);

        if ((time_counter % checkpoint_interval) == 0) {
            export_data(t, x, y, u0, num_zones, time_counter+1);
            printf("checkpoint%d.txt\n", time_counter / checkpoint_interval);
        }

        time_counter += 1;

        clock_t finish = clock();

        double mzps = (CLOCKS_PER_SEC / ((double)(finish - start))) *
                       num_zones * num_zones / 1E6;

        printf("[t=%.2e, dt=%.2e, Mzps=%.2e]\n", t, dt, mzps);
    }
    free(u1);
    free(u2);
    free(u3);
    free(sources);
}

int main()
{
    // initialize spatial grid, time, and concentration
    double x[num_zones];
    double y[num_zones];
    const double x_l = 0.0;
    const double x_r = 50.0;
    const double dx = (x_r - x_l) / num_zones;
    const double y_l = 0.0;
    const double y_r = 50.0;
    const double dy = (y_r - y_l) / num_zones;

    construct_grid(x, y, x_l, y_l, dx, dy);

    double t = 0.0;
    double *u0 = malloc(num_zones * num_zones * sizeof(double));
    double *ui = malloc(num_zones * num_zones * sizeof(double));
    for (int i = 0; i < num_zones; ++i) {
        for (int j = 0; j < num_zones; ++j) {
            u0[i * num_zones + j] = initial_condition(x[i], y[j]);
            ui[i * num_zones + j] = initial_condition(x[i], y[j]);
        }
    }

    // save t=0 checkpoint
    export_data(t, x, y, u0, num_zones, 0);

    // evolve the simulation in time
    rk3(u0, t, x, y, dx, dy);

    free(u0);
    free(ui);

    return 0;
}
