// -----------------------------------------------------------------------------
// Save data to text file.
// -----------------------------------------------------------------------------


// exports coordinates, initial u, and final u to a text file
void export_data(double *x, double *y, double *u, double *u0, int num_zones) {
  FILE *f = fopen("data/data.txt", "w");
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
