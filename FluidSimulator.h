#pragma once
#define ZERO_CELCIUS 273
#include "array2.h"
#include "pcgsolver/pcg_solver.h"
#include "pcgsolver/sparse_matrix.h"
#include "vec.h"
#include "MacGrid.h"
#include "LevelSet.h"

struct ScalarSource {
    float x, y, value;
};

struct VectorSource {
    float x, y;
    float val_x, val_y;
};

class FluidSimulator {
public:
    //Get velocity at point
    //Coordinates given on original grid

    //Get scalar quantity at point

    void advance(MacGrid& grid, float t_frame);
    void advance(MacGrid& grid, float t_frame, std::vector<ScalarSource>* temperatureSrcs, std::vector<ScalarSource>* concentrationSrcs, std::vector<VectorSource>* velocitySrcs, 
        bool isTemperatureSrcActive = false, bool isConcentrationSrcActive = false, bool isVelocitySrcActive = false);

    const float grav = 9.81f;
    const float density = 1;
    const float rate_t = 10;
    const float rate_c = 10;
    const float diffuse_cosnt_t = 1;
    const float ambient_temp = ZERO_CELCIUS + 21;
    const float diffuse_rate_temperature = 0.1;
    const float diffuse_rate_concentration = 0;
    FluidSimulator();
    ~FluidSimulator();
    void print_temperature_field(MacGrid& grid, const char* variable_name);
    void print_concentration_field(MacGrid& grid, const char* variable_name);
    void print_velocity_field(MacGrid& grid, const char* variable_name);
private:
    void get_advected_position(MacGrid& grid, float x_initial, float y_initial, float dt, float* x_result, float* y_result);

    float compute_cfl(MacGrid& grid);

    //void advect_particles(float dt);
    //Semi-Lagrangian advection
    void advect(MacGrid& grid, float dt);

    void advect_particles(MacGrid& grid, float dt);
    //Add external forces like gravity
    void add_forces(MacGrid& grid, float dt);

    //Apply incrompressibility constraint
    void project(MacGrid& grid, float dt);
    //void constrain_velocity();
    void enforceDirichlet(MacGrid& mac_grid);

    void classifyCells(MacGrid& grid, LevelSet& levelSet);
    void add_temperature_source(MacGrid& grid, const ScalarSource& src, float dt);
    void add_concentration_source(MacGrid& grid, const ScalarSource& src, float dt);
    void add_velocity_source(MacGrid& grid, const VectorSource& src, float dt);
    
    void diffuse_scalar(MacGrid& grid, float dt, float diffuse_rate);

    // Implement with std::function
    template <class T, class VectorT = std::vector<T>>
    void apply_sources(MacGrid& grid, VectorT* sources, float dt, void (FluidSimulator::* function)(MacGrid&, const T&, float));

    //void advectLevelSet(MacGrid& grid, LevelSet& levelSet, float dt);
    //helpers for pressure projection
    //void compute_weights();
    //void solve_pressure(float dt);
    //void compute_phi
};
