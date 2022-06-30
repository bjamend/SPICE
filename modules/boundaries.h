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
