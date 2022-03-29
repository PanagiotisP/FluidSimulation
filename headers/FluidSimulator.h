#pragma once
#include "FluidDomain.h"
#include "LevelSet.h"
#include "MacGrid.h"
#include "array2.h"
#include "pcgsolver/pcg_solver.h"
#include "pcgsolver/sparse_matrix.h"
#include "vec.h"

struct ScalarSource {
    float x, y, value;
};

struct VectorSource {
    float x, y;
    float val_x, val_y;
};

class FluidSimulator {
public:
    void transfer_from_particles_to_grid(FluidDomain &domain);
    void transfer_from_grid_to_particles(FluidDomain &domain, float flip_pic_ratio);

    void advance_flip_pic(FluidDomain &domain, float t_frame, float flip_pic_ratio);

    void advance_eulerian_grid(FluidDomain &domain, float t_frame);
    void advance_eulerian_grid(FluidDomain &domain, float t_frame, std::vector<ScalarSource> *temperatureSrcs,
                               std::vector<ScalarSource> *concentrationSrcs, std::vector<VectorSource> *velocitySrcs,
                               bool isTemperatureSrcActive = false, bool isConcentrationSrcActive = false,
                               bool isVelocitySrcActive = false);

    const float grav = 9.81f;

    FluidSimulator();
    ~FluidSimulator();
    void print_temperature_field(MacGrid &grid, const char *variable_name);
    void print_concentration_field(MacGrid &grid, const char *variable_name);
    void print_velocity_field(MacGrid &grid, const char *variable_name);

private:
    void get_advected_position(MacGrid &grid, float x_initial, float y_initial, float dt, float *x_result,
                               float *y_result);

    float compute_cfl(MacGrid &grid);

    // void advect_particles(float dt);
    // Semi-Lagrangian advection
    void advect(FluidDomain &domain, float dt);

    void advect_particles(FluidDomain &domain, float dt);
    void reseeding(FluidDomain &domain);
    // Add external forces like gravity
    void add_forces(FluidDomain &domain, float dt);

    // Apply incrompressibility constraint
    void project(FluidDomain &domain, float dt);
    // void constrain_velocity();
    void enforceDirichlet(FluidDomain &domain);

    void add_temperature_source(MacGrid &grid, const ScalarSource &src, float dt);
    void add_concentration_source(MacGrid &grid, const ScalarSource &src, float dt);
    void add_velocity_source(MacGrid &grid, const VectorSource &src, float dt);

    void diffuse_scalar(FluidDomain &domain, float dt, float diffuse_rate);

    // Implement with std::function
    template <class T, class VectorT = std::vector<T>>
    void apply_sources(MacGrid &grid, VectorT *sources, float dt,
                       void (FluidSimulator::*function)(MacGrid &, const T &, float));
    
    float calculate_kernel_function(MacGrid& grid, float x, float y);

    void extrapolate_data(FluidDomain& domain, int iterations_n);
};
