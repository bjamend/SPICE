// -----------------------------------------------------------------------------
// Input parameters for simulation
// -----------------------------------------------------------------------------


const int num_zones = 256;   // number of spatial zones along one axis
const double t_final = 0.01; // time at which to end the simulation [Gyr]
const double cfl = 0.4;      // value to satisfy CFL convergence condition
const double r = 0.01;       // initial radius of r-process event [kpc]
const double source_production_rate = 1e5;
const int max_sources_per_timestep = 1000;
const double source_amplitude = 10000.0;
