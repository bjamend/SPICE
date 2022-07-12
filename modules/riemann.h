// -----------------------------------------------------------------------------
// Compute fluxes to solve the Riemann problem (2nd-order).
// -----------------------------------------------------------------------------


#include "model.h"
#include "extras.h"
#include "grid.h"


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
double source_terms(double x, double y, double t, double *events) {
  double result = 0.0;
  for (int i = 0; i < 10; ++i) {
    double t_0 = events[3*i];
    double x_0 = events[3*i+1];
    double y_0 = events[3*i+2];
    result += (exp(-((x - x_0) * (x - x_0) + (y - y_0) * (y - y_0)) / r / r) *
               exp(-(t - t_0) * (t - t_0) / tau / tau));
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
             double t, double *events) {

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
  double source = source_terms(x, y, t, events);
  return -(f_iphj - f_imhj) / dx - (g_ijph - g_ijmh) / dy + source;
}
