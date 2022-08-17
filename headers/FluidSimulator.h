#pragma once
#include "FluidDomain.h"
#include "LevelSet.h"
#include "MacGrid.h"
#include "SolidObject.h"
#include "pcgsolver/pcg_solver.h"
#include "pcgsolver/sparse_matrix.h"
#include "vec.h"

#include <openvdb/openvdb.h>

class FluidSimulator {
public:
    FluidSimulator();
    FluidSimulator(openvdb::CoordBBox bbox);
    ~FluidSimulator();

    void advance_flip_pic(FluidDomain &domain, float t_frame, float flip_pic_ratio);
    const float grav = 9.81f;

    void print_velocity_field(MacGrid &grid, const char *variable_name);

    static double calculate_kernel_function(openvdb::Vec3d vec);
    static double calculate_kernel_function(double x, double y, double z);
    static openvdb::Vec3d calculate_kernel_function_staggered(openvdb::Vec3d difference);

private:
    float compute_cfl(FluidDomain &domain);
    void transfer_from_particles_to_grid(FluidDomain &domain);
    void reseeding(FluidDomain &domain);
    // Add external forces like gravity
    void add_forces(FluidDomain &domain, float dt);
    void extrapolate_data(FluidDomain &domain, int iterations_n);
    void constrain_velocity(FluidDomain &domain);
    void constrain_velocity(FluidDomain &domain, SolidObject solidObject);
    void project(FluidDomain &domain, float dt);
    // Apply incrompressibility constraint
    // Semi-Lagrangian advection
    void advect(FluidDomain &domain, float dt);
    void advect_particles(FluidDomain &domain, float dt);

    void transfer_from_grid_to_particles(FluidDomain &domain, float flip_pic_ratio);

    const openvdb::CoordBBox _bbox;
};
