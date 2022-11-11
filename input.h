// -----------------------------------------------------------------------------
// Input parameters for simulation
// -----------------------------------------------------------------------------


const int num_zones = 128;        // number of spatial zones along one axis
const double t_final = 10.0;      // time at which to end the simulation [Gyr]
const double cfl = 0.4;           // value to satisfy CFL convergence condition
const double r = 0.1;             // initial radius of r-process event [kpc]
const double source_production_rate = 1e4;
const int max_sources_per_timestep = 1e8;
const double source_yield = 1.06E3;  // yield produced by single source [M.]
const int checkpoint_interval = 5E0;
