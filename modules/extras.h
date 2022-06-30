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
