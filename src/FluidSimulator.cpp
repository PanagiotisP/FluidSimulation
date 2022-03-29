#pragma once
#include "FluidSimulator.h"

#include "FluidDomain.h"
#include "array2_utils.h"
#include "util.h"
#include "vec.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <random>

// #define PRINT
#define VERTICAL_PLANE true
const float epsilon = 10e-37;

float liquid_bound_x0, liquid_bound_x1;
float liquid_bound_y0, liquid_bound_y1;

float solid_bound_x0, solid_bound_x1;
float solid_bound_y0, solid_bound_y1;


float distance_from_axis_aligned_rect_box(const float x0, const float x1, const float y0, const float y1,
                                          const Vec2f &point) {
    // inside
    if (x0 < point[0] && point[0] < x1 && y0 < point[1] && point[1] < y1) {
        return max(x0 - point[0], point[0] - x1, y0 - point[1], point[1] - y1);
    }
    // outside
    else {
        // closest point to bounding box (p,q)
        float p, q;
        if (point[0] < x0)
            p = x0;
        else if (x1 < point[0])
            p = x1;
        else
            p = point[0];

        if (point[1] < y0)
            q = y0;
        else if (y1 < point[1])
            q = y1;
        else
            q = point[1];
        return sqrt((sqr(point[0] - p)) + (sqr(point[1] - q)));
    }
}

FluidSimulator::FluidSimulator() {}

FluidSimulator::~FluidSimulator() {}

void FluidSimulator::print_temperature_field(MacGrid &grid, const char *variable_name) {
#ifdef PRINT
    std::cout << variable_name;
    for (int j = grid.sizeY() - 1; j >= 0; --j) {
        std::cout << std::endl;
        for (int i = 0; i < grid.sizeX(); ++i) { std::cout << grid.temperature(i, j) - ZERO_CELCIUS << " "; }
    }
    std::cout << std::endl;
#endif // PRINT
}

void FluidSimulator::print_concentration_field(MacGrid &grid, const char *variable_name) {
#ifdef PRINT
    std::cout << variable_name;
    for (int j = grid.sizeY() - 1; j >= 0; --j) {
        std::cout << std::endl;
        for (int i = 0; i < grid.sizeX(); ++i) { std::cout << grid.concentration(i, j) << " "; }
    }
    std::cout << std::endl;
#endif // PRINT
}

void FluidSimulator::print_velocity_field(MacGrid &grid, const char *variable_name) {
#ifdef PRINT
    std::cout << variable_name;
    for (int j = grid.sizeY() - 1; j >= 0; --j) {
        std::cout << std::endl;
        for (int i = 0; i < grid.sizeX(); ++i) {
            std::cout << '[' << grid.velXHalfIndexed(i, j) << ", " << grid.velYHalfIndexed(i, j) << "] ";
        }
    }
    std::cout << std::endl;
#endif // PRINT
}

void FluidSimulator::diffuse_scalar(FluidDomain &domain, float dt, float diffuse_rate) {
    auto &grid = domain.grid();
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (domain.fluidLevelSet()(i, j) < 0) {
                int i_minus1 = clamp(i - 1, 0, grid.sizeX() - 1);
                int j_minus1 = clamp(j - 1, 0, grid.sizeY() - 1);

                int i_plus1 = clamp(i + 1, 0, grid.sizeX() - 1);
                int j_plus1 = clamp(j + 1, 0, grid.sizeY() - 1);


                float neighbouring_values = (grid.concentration(i_minus1, j) + grid.concentration(i, j_minus1)
                                             + grid.concentration(i_plus1, j) + grid.concentration(i, j_plus1))
                                            / 4;
                float next_value = (grid.concentration(i, j) + diffuse_rate * neighbouring_values) / (1 + diffuse_rate);

                grid.setConcentrationBackBuffer(i, j, next_value);
            }
        }
    }
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) { grid.setConcentration(i, j, grid.concentrationBackBuffer(i, j)); }
    }
}

void FluidSimulator::add_temperature_source(MacGrid &grid, const ScalarSource &src, float dt) {
    float current_temperature = grid.temperatureInterpolated(src.x, src.y);
    float temperature_increase = (1 - exp(-FluidDomain::rate_t * dt)) * (src.value - current_temperature);
    grid.addToTemperatureInterpolated(src.x, src.y, temperature_increase);
}

void FluidSimulator::add_concentration_source(MacGrid &grid, const ScalarSource &src, float dt) {
    float current_concentration = grid.concentrationInterpolated(src.x, src.y);
    float concentration_increase = dt * (src.value - current_concentration);
    grid.addToConcentrationInterpolated(src.x, src.y, concentration_increase);
}

void FluidSimulator::add_velocity_source(MacGrid &grid, const VectorSource &src, float dt) {
    grid.addToVelXInterpolated(src.x, src.y, src.val_x);
    grid.addToVelYInterpolated(src.x, src.y, src.val_y);
}

template <class T, class VectorT>
void FluidSimulator::apply_sources(MacGrid &grid, VectorT *sources, float dt,
                                   void (FluidSimulator::*function)(MacGrid &, const T &, float)) {
    for (auto it = sources->begin(); it != sources->end(); ++it) { (this->*function)(grid, *it, dt); }
}

void FluidSimulator::transfer_from_particles_to_grid(FluidDomain &domain) {
    MacGrid &grid = domain.grid();
    Grid<float> vel_x_weight(grid.sizeX(), grid.sizeY(), grid.deltaX(), grid.deltaY());
    Grid<float> vel_y_weight(grid.sizeX(), grid.sizeY(), grid.deltaX(), grid.deltaY());
    Grid<float> center_weight(grid.sizeX(), grid.sizeY(), grid.deltaX(), grid.deltaY());

    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            grid.setVelXBackBufferHalfIndexed(i, j, 0);
            grid.setVelYBackBufferHalfIndexed(i, j, 0);
            grid.setTemperatureBackBuffer(i, j, FluidDomain::ambient_temp);
            grid.setConcentrationBackBuffer(i, j, 0);

            vel_x_weight(i, j) = 0;
            vel_y_weight(i, j) = 0;
            center_weight(i, j) = 0;
        }
    }
    float k;

    // The domain of kernel function that gives non-zero value is [-1.5, 1.5]
    // That means that for a given cell (i,j) the grid points with a possible non-zero value are 16:
    // The corners of the grid cell + the corners of the adjacent in 8-way connectivity grid cells
    //  o -   o   - o - o
    //  o -   o   - o - o
    //  o - (i,j) - o - o
    //  o -   o   - o - o
    for (auto it = domain.particleSet().begin(); it != domain.particleSet().end(); ++it) {
        auto &p = *it;
        int i_start = clamp(static_cast<int>(floor(p->posX() / domain.grid().deltaX())) - 1, 0, grid.sizeX() - 1);
        int j_start = clamp(static_cast<int>(floor(p->posY() / domain.grid().deltaY())) - 1, 0, grid.sizeY() - 1);

        int i_end = clamp(i_start + 4, 0, grid.sizeX() - 1);
        int j_end = clamp(j_start + 4, 0, grid.sizeY() - 1);


        for (int j = j_start; j < j_end; ++j) {
            for (int i = i_start; i < i_end; ++i) {
                // velocity x-component
                k = calculate_kernel_function(grid, p->posX() - i * grid.deltaX(),
                                              p->posY() - (j + 0.5) * grid.deltaY());
                vel_x_weight(i, j) += k;
                grid.setVelXBackBufferHalfIndexed(i, j, grid.velXBackBufferHalfIndexed(i, j) + k * p->velX());

                // velocity y-component
                k = calculate_kernel_function(grid, p->posX() - (i + 0.5) * grid.deltaX(),
                                              p->posY() - j * grid.deltaY());
                vel_y_weight(i, j) += k;
                grid.setVelYBackBufferHalfIndexed(i, j, grid.velYBackBufferHalfIndexed(i, j) + (k * p->velY()));


                // quantities sampled at the center of the cell
                k = calculate_kernel_function(grid, p->posX() - (i + 0.5) * grid.deltaX(),
                                              p->posY() - (j + 0.5) * grid.deltaY());
                center_weight(i, j) += k;
                // temperature
                grid.setTemperatureBackBuffer(i, j, grid.temperatureBackBuffer(i, j) + (k * p->temperature()));

                // concentration
                grid.setConcentrationBackBuffer(i, j, grid.concentrationBackBuffer(i, j) + (k * p->concentration()));
            }
        }
    }
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (vel_x_weight(i, j) > 0) {
                grid.setVelXBackBufferHalfIndexed(i, j, grid.velXBackBufferHalfIndexed(i, j) / vel_x_weight(i, j));
            }
            if (vel_y_weight(i, j) > 0) {
                grid.setVelYBackBufferHalfIndexed(i, j, grid.velYBackBufferHalfIndexed(i, j) / vel_y_weight(i, j));
            }
            if (center_weight(i, j) > 0) {
                grid.setTemperatureBackBuffer(i, j, grid.temperatureBackBuffer(i, j) / center_weight(i, j));
                grid.setConcentrationBackBuffer(i, j, grid.concentrationBackBuffer(i, j) / center_weight(i, j));
            }
        }
    }
    grid.swapVelocityBuffers();
    grid.swapTemperatureBuffers();
    grid.swapConcentrationBuffers();
}

void FluidSimulator::transfer_from_grid_to_particles(FluidDomain &domain, float flip_pic_ratio = 0.98) {
    MacGrid &grid = domain.grid();
    for (auto it = domain.particleSet().begin(); it != domain.particleSet().end(); ++it) {
        auto &p = *it;
        float pic_vel_x = grid.velXInterpolated(p->posX(), p->posY());
        float pic_vel_y = grid.velYInterpolated(p->posX(), p->posY());
        float pic_temperature = grid.temperatureInterpolated(p->posX(), p->posY());
        float pic_concentration = grid.concentrationInterpolated(p->posX(), p->posY());

        float flip_vel_x = p->velX() + grid.velXDiffInterpolated(p->posX(), p->posY());
        float flip_vel_y = p->velY() + grid.velYDiffInterpolated(p->posX(), p->posY());
        float flip_temperature = p->temperature() + grid.temperatureDiffInterpolated(p->posX(), p->posY());
        float flip_concentration = p->concentration() + grid.concentrationDiffInterpolated(p->posX(), p->posY());


        p->setVelocity(pic_vel_x * (1 - flip_pic_ratio) + flip_vel_x * flip_pic_ratio,
                       pic_vel_y * (1 - flip_pic_ratio) + flip_vel_y * flip_pic_ratio);
        p->setTemperature(pic_temperature * (1 - flip_pic_ratio) + flip_temperature * flip_pic_ratio);
        p->setConcentration(pic_concentration * (1 - flip_pic_ratio) + flip_concentration * flip_pic_ratio);
    }
}

void FluidSimulator::reseeding(FluidDomain &domain) {
    MacGrid &grid(domain.grid());
    std::vector<ParticleSet::iterator> toBeErased;
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<> dis(0, 1);
    auto compareParticles = [](const ParticleSet::iterator &p1, const ParticleSet::iterator &p2) {
        return (*p1)->concentration() < (*p2)->concentration();
    };
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (domain.fluidLevelSet()(i, j) < 0) {
                int particle_counter = 0;
                int k = 0;
                toBeErased.clear();
                for (auto it = domain.particleSet().begin(); it != domain.particleSet().end(); ++it) {
                    auto &p = *it;
                    if (i * grid.deltaX() < p->posX() && p->posX() <= (i + 1) * grid.deltaX()
                        && j * grid.deltaY() < p->posY() && p->posY() <= (j + 1) * grid.deltaY()) {
                        particle_counter++;
                        toBeErased.push_back(it);
                    }
                }
                if (particle_counter > 8) {
                    // Remove the particles based on concentration
                    std::sort(toBeErased.begin(), toBeErased.end(), compareParticles);
                    toBeErased.resize(particle_counter - 4);
                    std::sort(toBeErased.begin(), toBeErased.end());
                    while (toBeErased.size() > 8) {
                        domain.particleSet().removeParticle(toBeErased.back());
                        toBeErased.pop_back();
                    }
                } else if (particle_counter == 1) {
                    float new_x = (i + dis(gen)) * grid.deltaX();
                    float new_y = (j + dis(gen)) * grid.deltaY();

                    int max_tries = 100;
                    while (domain.fluidLevelSet().valueInterpolated(new_x - 0.5 * grid.deltaX(),
                                                                    new_y - 0.5 * grid.deltaY())
                               > 0
                           && max_tries-- >= 0) {
                        new_x = (i + dis(gen)) * grid.deltaX();
                        new_y = (j + dis(gen)) * grid.deltaY();
                    }
                    // Limit exceeded
                    if (max_tries >= 0) {
                        auto p = std::make_unique<Particle>(
                            new_x, new_y, grid.velXInterpolated(new_x, new_y), grid.velYInterpolated(new_x, new_y),
                            grid.temperatureInterpolated(new_x, new_y), grid.concentrationInterpolated(new_x, new_y));
                        domain.particleSet().addParticle(p);
                    }
                }
            }
        }
    }
}

void FluidSimulator::extrapolate_data(FluidDomain &domain, int iterations_n) {
    auto &grid = domain.grid();
    Grid<int> valid_mask_x_front_buffer(grid.sizeX(), grid.sizeY(), grid.deltaX(), grid.deltaY());
    Grid<int> valid_mask_x_back_buffer(grid.sizeX(), grid.sizeY(), grid.deltaX(), grid.deltaY());
    Grid<int> valid_mask_y_front_buffer(grid.sizeX(), grid.sizeY(), grid.deltaX(), grid.deltaY());
    Grid<int> valid_mask_y_back_buffer(grid.sizeX(), grid.sizeY(), grid.deltaX(), grid.deltaY());

    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            int i_minus_1 = clamp(i - 1, 0, grid.sizeX() - 1);
            int j_minus_1 = clamp(j - 1, 0, grid.sizeY() - 1);
            if (domain.fluidLevelSet()(i, j) < 0 || domain.fluidLevelSet()(i_minus_1, j) < 0) {
                valid_mask_x_front_buffer(i, j) = 1;
                valid_mask_x_back_buffer(i, j) = 1;
                grid.setVelXBackBufferHalfIndexed(i, j, grid.velXHalfIndexed(i, j));
            } else {
                valid_mask_x_front_buffer(i, j) = 0;
                valid_mask_x_back_buffer(i, j) = 0;
                grid.setVelXBackBufferHalfIndexed(i, j, 0);
                grid.setVelXHalfIndexed(i, j, 0);
            }
            if (domain.fluidLevelSet()(i, j) < 0 || domain.fluidLevelSet()(i, j_minus_1) < 0) {
                valid_mask_y_front_buffer(i, j) = 1;
                valid_mask_y_back_buffer(i, j) = 1;
                grid.setVelYBackBufferHalfIndexed(i, j, grid.velYHalfIndexed(i, j));
            } else {
                valid_mask_y_front_buffer(i, j) = 0;
                valid_mask_y_back_buffer(i, j) = 0;
                grid.setVelYBackBufferHalfIndexed(i, j, 0);
                grid.setVelYHalfIndexed(i, j, 0);
            }
        }
    }

    for (int iter = 0; iter < iterations_n; ++iter) {
        for (int j = 0; j < grid.sizeY(); ++j) {
            for (int i = 0; i < grid.sizeX(); ++i) {
                int i_minus_1 = clamp(i - 1, 0, grid.sizeX() - 1);
                int j_minus_1 = clamp(j - 1, 0, grid.sizeY() - 1);
                int i_plus_1 = clamp(i + 1, 0, grid.sizeX() - 1);
                int j_plus_1 = clamp(j + 1, 0, grid.sizeY() - 1);
                // Check if this cell will be updated in the x dimension
                if (valid_mask_x_front_buffer(i, j) == 0 && domain.solidLevelSet()(i, j) > 0
                    && domain.solidLevelSet()(i_minus_1, j) > 0) {
                    float new_vel_x = 0;
                    int n_valid_neighbors_x = 0;

                    // Get values of all valid neighbors
                    if (valid_mask_x_front_buffer(i_minus_1, j) == 1) {
                        new_vel_x += grid.velXBackBufferHalfIndexed(i_minus_1, j);
                        n_valid_neighbors_x++;
                    }
                    if (valid_mask_x_front_buffer(i, j_minus_1) == 1) {
                        new_vel_x += grid.velXBackBufferHalfIndexed(i, j_minus_1);
                        n_valid_neighbors_x++;
                    }
                    if (valid_mask_x_front_buffer(i, j_plus_1) == 1) {
                        new_vel_x += grid.velXBackBufferHalfIndexed(i, j_plus_1);
                        n_valid_neighbors_x++;
                    }
                    if (valid_mask_x_front_buffer(i_plus_1, j) == 1) {
                        new_vel_x += grid.velXBackBufferHalfIndexed(i_plus_1, j);
                        n_valid_neighbors_x++;
                    }

                    // Average the value for the current cell
                    if (n_valid_neighbors_x > 0) {
                        new_vel_x /= n_valid_neighbors_x;
                        grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);
                        valid_mask_x_back_buffer(i, j) = 1;
                    }
                }

                // Check if this cell will be updated in the y dimension
                if (valid_mask_y_front_buffer(i, j) == 0 && domain.solidLevelSet()(i, j) > 0
                    && domain.solidLevelSet()(i, j_minus_1) > 0) {
                    float new_vel_y = 0;
                    int n_valid_neighbors_y = 0;

                    // Get values of all valid neighbors
                    if (valid_mask_y_front_buffer(i_minus_1, j) == 1) {
                        new_vel_y += grid.velYBackBufferHalfIndexed(i_minus_1, j);
                        n_valid_neighbors_y++;
                    }
                    if (valid_mask_y_front_buffer(i, j_minus_1) == 1) {
                        new_vel_y += grid.velYBackBufferHalfIndexed(i, j_minus_1);
                        n_valid_neighbors_y++;
                    }
                    if (valid_mask_y_front_buffer(i, j_plus_1) == 1) {
                        new_vel_y += grid.velYBackBufferHalfIndexed(i, j_plus_1);
                        n_valid_neighbors_y++;
                    }
                    if (valid_mask_y_front_buffer(i_plus_1, j) == 1) {
                        new_vel_y += grid.velYBackBufferHalfIndexed(i_plus_1, j);
                        n_valid_neighbors_y++;
                    }

                    // Average the value for the current cell
                    if (n_valid_neighbors_y > 0) {
                        new_vel_y /= n_valid_neighbors_y;
                        grid.setVelYBackBufferHalfIndexed(i, j, new_vel_y);
                        valid_mask_y_back_buffer(i, j) = 1;
                    }
                }
            }
        }
        // Swap valid mask
        valid_mask_x_front_buffer = valid_mask_x_back_buffer;
        valid_mask_y_front_buffer = valid_mask_y_back_buffer;
    }
    grid.swapVelocityBuffers();
}

void FluidSimulator::advance_flip_pic(FluidDomain &domain, float t_frame, float flip_pic_ratio = 0.98) {
    float t = 0;
    while (t < t_frame) {
        float timestep = compute_cfl(domain.grid());
        if (timestep + t >= t_frame) { timestep = t_frame - t; }
        printf("Taking substep of size %f (to %0.3f%% of the frame)\n", timestep, 100 * (t + timestep) / t_frame);
        domain.update(timestep);
        domain.construct_level_set_from_marker_particles(domain.fluidLevelSet());

        // reseeding(domain);

        transfer_from_particles_to_grid(domain);
        domain.grid().updatePreviousVelocityBuffer();

        // Eulerian grid part START
        add_forces(domain, timestep);
        enforceDirichlet(domain);
        extrapolate_data(domain, 2);
        project(domain, timestep);
        enforceDirichlet(domain);
        // Eulerian grid part END

        domain.grid().updateDiffBuffers();

        transfer_from_grid_to_particles(domain, flip_pic_ratio);
        domain.advectParticles(timestep);

        t += timestep;
    }
}

void FluidSimulator::advance_eulerian_grid(FluidDomain &domain, float t_frame) {
    std::vector<ScalarSource> *temp_scalar_vec = nullptr;
    std::vector<VectorSource> *temp_vector_vec = nullptr;
    advance_eulerian_grid(domain, t_frame, temp_scalar_vec, temp_scalar_vec, temp_vector_vec, false, false, false);
}

void FluidSimulator::advance_eulerian_grid(FluidDomain &domain, float t_frame,
                                           std::vector<ScalarSource> *temperatureSrcs,
                                           std::vector<ScalarSource> *concentrationSrcs,
                                           std::vector<VectorSource> *velocitySrcs, bool isTemperatureSrcActive,
                                           bool isConcentrationSrcActive, bool isVelocitySrcActive) {
    float t = 0;
    while (t < t_frame) {
        float timestep = compute_cfl(domain.grid());
        if (timestep + t >= t_frame) { timestep = t_frame - t; }
        printf("Taking substep of size %f (to %0.3f%% of the frame)\n", timestep, 100 * (t + timestep) / t_frame);

        advect(domain, timestep);
        print_velocity_field(domain.grid(), "");

        if (isTemperatureSrcActive) {
            apply_sources<ScalarSource>(domain.grid(), temperatureSrcs, timestep,
                                        &FluidSimulator::add_temperature_source);
            domain.grid().swapTemperatureBuffers();
        }

        if (isConcentrationSrcActive) {
            apply_sources<ScalarSource>(domain.grid(), concentrationSrcs, timestep,
                                        &FluidSimulator::add_concentration_source);
            domain.grid().swapConcentrationBuffers();
        }

        if (isVelocitySrcActive) {
            apply_sources<VectorSource>(domain.grid(), velocitySrcs, timestep, &FluidSimulator::add_velocity_source);
            domain.grid().swapVelocityBuffers();
        }
        print_velocity_field(domain.grid(), "");


        diffuse_scalar(domain, timestep, 0.01);

        // apply forces (e.g. gravity)
        add_forces(domain, timestep);
        // print_velocity_field(domain.grid(), "After forces");
        // print_temperature_field(domain.grid(), "After forces");


        project(domain, timestep);
        print_velocity_field(domain.grid(), "");
        // print_velocity_field(domain.grid(), "After projection");

        enforceDirichlet(domain);

        ////update temperature and concentration
        // float rate_t = 10;
        // float rate_s = 1;
        // float target_temperature = 273 + 37;
        ////temperature(domain.grid().sizeX() / 2, domain.grid().sizeY() / 2) += (1 - exp(-rate_t * timestep)) *
        ///(target_temperature - temperature(domain.grid().sizeX() / 2, domain.grid().sizeY() / 2));
        // temperature(0, 0) += (1 - exp(-rate_t * timestep)) * (target_temperature - temperature(0, 0));
        // temperature(0, 1) += (1 - exp(-rate_t * timestep)) * (target_temperature - temperature(0, 1));
        // concentration(domain.grid().sizeX() / 2, domain.grid().sizeY() / 2) += clamp(rate_s * timestep,
        // 0.f, 1.f);

        t += timestep;
    }
}

void FluidSimulator::get_advected_position(MacGrid &grid, float x_initial, float y_initial, float dt, float *x_result,
                                           float *y_result) {
    float x_mid = x_initial - 0.5 * dt * grid.velXInterpolated(x_initial, y_initial);
    float y_mid = y_initial - 0.5 * dt * grid.velYInterpolated(x_initial, y_initial);

    *x_result = x_initial - dt * grid.velXInterpolated(x_mid, y_mid);
    *y_result = y_initial - dt * grid.velYInterpolated(x_mid, y_mid);
    return;
}

// void FluidSimulator::advectLevelSet(MacGrid& grid, LevelSet& levelSet, float dt) {
//     LevelSet newlevelSet(
//         levelSet.sizeX(),
//         levelSet.sizeY(),
//         levelSet.lengthX(),
//         levelSet.lengthY());
//     for (int j = 0; j < levelSet.sizeX(); ++j) 	{
//         for (int i = 0; i < levelSet.sizeY(); ++i) 		{
//             float vel_x = grid.velX(i, j);
//             float vel_y = grid.velY(i, j);
//             float grad_x = levelSet.computeUpwindGradientX(i, j, vel_x);
//             float grad_y = levelSet.computeUpwindGradientY(i, j, vel_y);
//
//             float change_rate = -(vel_x * grad_x + vel_y * grad_y);
//
//             // Forward Euler
//             newlevelSet(i, j) = levelSet(i, j) + change_rate * dt * 5;
//         }
//     }
//     for (int j = 0; j < levelSet.sizeX(); ++j) {
//         for (int i = 0; i < levelSet.sizeY(); ++i) {
//             levelSet(i, j) = newlevelSet(i, j);
//         }
//     }
// }

float FluidSimulator::compute_cfl(MacGrid &grid) {
    double max_vel = 0;
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            max_vel = max((double)fabs(grid.velXHalfIndexed(i, j)), max_vel);
            max_vel = max((double)fabs(grid.velYHalfIndexed(i, j)), max_vel);
        }
    }
    max_vel += sqrt(5 * max(grid.deltaX(), grid.deltaY()) * grav);
    return min(grid.deltaX(), grid.deltaY()) / (max_vel + epsilon);
}

void FluidSimulator::advect(FluidDomain &domain, float dt) {
    auto &grid = domain.grid();
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            grid.setVelXBackBufferHalfIndexed(i, j, grid.velXHalfIndexed(i, j));
            grid.setVelYBackBufferHalfIndexed(i, j, grid.velYHalfIndexed(i, j));
            grid.setTemperatureBackBuffer(i, j, FluidDomain::ambient_temp);
            grid.setConcentrationBackBuffer(i, j, 0);
        }
    }

    // Hypothetical particle for x-component
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (domain.fluidLevelSet()(i, j) < 0 || domain.fluidLevelSet()(i - 1, j) < 0) {
                float x_current = i * grid.deltaX();
                float y_current = (j + 0.5) * grid.deltaY();
                float x_previous, y_previous;

                get_advected_position(grid, x_current, y_current, dt, &x_previous, &y_previous);

                float v_x = grid.velXInterpolated(x_previous, y_previous);
                grid.setVelXBackBufferHalfIndexed(i, j, v_x);
            }
        }
    }

    // Hypothetical particle for y-component
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (domain.fluidLevelSet()(i, j) < 0 || domain.fluidLevelSet()(i - 1, j) < 0) {
                float x_current = (i + 0.5) * grid.deltaX();
                float y_current = j * grid.deltaY();
                float x_previous, y_previous;

                get_advected_position(grid, x_current, y_current, dt, &x_previous, &y_previous);

                float v_y = grid.velYInterpolated(x_previous, y_previous);
                grid.setVelYBackBufferHalfIndexed(i, j, v_y);
            }
        }
    }

    // Advect other quantities like Temperature and Concentration
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (domain.fluidLevelSet()(i, j) < 0 || domain.fluidLevelSet()(i - 1, j) < 0) {
                float x_current = (i + 0.5) * grid.deltaX();
                float y_current = (j + 0.5) * grid.deltaY();
                float x_previous, y_previous;

                get_advected_position(grid, x_current, y_current, dt, &x_previous, &y_previous);

                float temperature = grid.temperatureInterpolated(x_previous, y_previous);
                grid.setTemperatureBackBuffer(i, j, temperature);
                float concentration = grid.concentrationInterpolated(x_previous, y_previous);
                grid.setConcentrationBackBuffer(i, j, concentration);
            }
        }
    }

    // Write results to front buffer
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            grid.setVelXHalfIndexed(i, j, grid.velXBackBufferHalfIndexed(i, j));
            grid.setVelYHalfIndexed(i, j, grid.velYBackBufferHalfIndexed(i, j));

            // Advect scalar field too
            grid.setTemperature(i, j, grid.temperatureBackBuffer(i, j));
            grid.setConcentration(i, j, grid.concentrationBackBuffer(i, j));
        }
    }
}

void FluidSimulator::add_forces(FluidDomain &domain, float dt) {
    // gravity
    auto &grid = domain.grid();
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            // Only add force to the liquid cells
            if (domain.fluidLevelSet()(i, j) < 0) {
                if (VERTICAL_PLANE) { grid.setVelYHalfIndexed(i, j, grid.velYHalfIndexed(i, j) - grav * dt); }
            }
        }
    }
}

void FluidSimulator::project(FluidDomain &domain, float dt) {
    auto &grid = domain.grid();
    int system_size = 0;

    // Solver of linear system
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float, 1 > > solver;
    // Laplacian matrix
    Eigen::SparseMatrix<float, 1> alpha_matrix;

    Grid<int> fluid_indices(grid.sizeX(), grid.sizeX(), grid.deltaX(), grid.deltaY());
    // Index all cells with fluid
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (domain.fluidLevelSet()(i, j) < 0) {
                fluid_indices(i, j) = system_size;
                system_size++;
            } else {
                fluid_indices(i, j) = -1;
            }
        }
    }
    if (system_size == 0) { return; }
    alpha_matrix.resize(system_size, system_size);
    alpha_matrix.reserve(Eigen::VectorXi::Constant(system_size, 5));

    Eigen::VectorXf rhs(system_size);

    // Pressure solver
    // std::vector<double> rhs(system_size);
    // std::vector<double> pressure_grid(system_size);
    // SparseMatrixd alpha_matrix(system_size);
    // PCGSolver<double> solver;

    // alpha_matrix.zero();
    // rhs.assign(rhs.size(), 0);
    // pressure_grid.assign(pressure_grid.size(), 0);


    float scale_x = dt / (FluidDomain::density * sqr(grid.deltaX()));
    float scale_y = dt / (FluidDomain::density * sqr(grid.deltaY()));

    // solid velocity zero for now
    float u_solid = 0;
    float v_solid = 0;

    // 2-d Array. Don't really need a grid for that
    Array2f u_weights(grid.sizeX(), grid.sizeY());
    Array2f v_weights(grid.sizeX(), grid.sizeY());

    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            float j_plus_1 = clamp(j + 1, 0, grid.sizeY() - 1);

            // TODO NEED TO REEVAULATE THE 0.5 (right now taking bottom left corners)
            float phi_i_j =
                domain.solidLevelSet().valueInterpolated((i - 0.5f) * grid.deltaX(), (j - 0.5f) * grid.deltaY());
            float phi_i_j_plus_1 =
                domain.solidLevelSet().valueInterpolated((i - 0.5f) * grid.deltaX(), (j_plus_1 - 0.5f) * grid.deltaY());

            sort(phi_i_j, phi_i_j_plus_1);
            if (phi_i_j == 0 && phi_i_j_plus_1 == 0) {
                u_weights(i, j) = 0;
            } else if (phi_i_j * phi_i_j_plus_1 <= 0) { // Solid - Fluid Boundary
                u_weights(i, j) = -phi_i_j_plus_1 / (phi_i_j - phi_i_j_plus_1);
            } else if (phi_i_j < 0) { // Completely in solid
                u_weights(i, j) = 0;
            } else { // Completely in fluid
                u_weights(i, j) = 1;
            }
            if (u_weights(i, j) < 0.1) u_weights(i, j) = 0;
        }
    }
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            float i_plus_1 = clamp(i + 1, 0, grid.sizeX() - 1);

            // TODO NEED TO REEVAULATE THE 0.5 (right now taking bottom left corners)
            float phi_i_j =
                domain.solidLevelSet().valueInterpolated((i - 0.5f) * grid.deltaX(), (j - 0.5f) * grid.deltaY());
            float phi_i_plus_1_j =
                domain.solidLevelSet().valueInterpolated((i_plus_1 - 0.5f) * grid.deltaX(), (j - 0.5f) * grid.deltaY());

            sort(phi_i_j, phi_i_plus_1_j);
            if (phi_i_j == 0 && phi_i_plus_1_j == 0) {
                v_weights(i, j) = 0;
            } else if (phi_i_j * phi_i_plus_1_j <= 0) { // Solid - Fluid Boundary
                v_weights(i, j) = -phi_i_plus_1_j / (phi_i_j - phi_i_plus_1_j);
            } else if (phi_i_j < 0) { // Completely in solid
                v_weights(i, j) = 0;
            } else { // Completely in fluid
                v_weights(i, j) = 1;
            }
            if (v_weights(i, j) < 0.1) v_weights(i, j) = 0;
        }
    }
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (fluid_indices(i, j) != -1) {
                int idx = fluid_indices(i, j); // i + j * grid.sizeX();
                rhs[idx] = 0;
                float diagonal = 0;
                float theta = 0;

                int i_minus_1 = clamp(i - 1, 0, grid.sizeX() - 1);
                int j_minus_1 = clamp(j - 1, 0, grid.sizeY() - 1);

                int i_plus_1 = clamp(i + 1, 0, grid.sizeX() - 1);
                int j_plus_1 = clamp(j + 1, 0, grid.sizeY() - 1);

                float phi_i_j = domain.fluidLevelSet()(i, j);
                float phi_i_minus_1_j = domain.fluidLevelSet()(i_minus_1, j);
                float phi_i_j_minus_1 = domain.fluidLevelSet()(i, j_minus_1);
                float phi_i_plus_1_j = domain.fluidLevelSet()(i_plus_1, j);
                float phi_i_j_plus_1 = domain.fluidLevelSet()(i, j_plus_1);

                // left neighbour
                // if (grid.cellType(i - 1, j) != SOLID) {
                if (i - 1 >= 0) {
                    if (phi_i_minus_1_j < 0) {
                        alpha_matrix.insert(idx, fluid_indices(i_minus_1, j)) = -u_weights(i, j) * scale_x;
                    }
                    theta = max(0.01f, fraction_inside(phi_i_minus_1_j, phi_i_j));
                    diagonal += u_weights(i, j) * scale_x * (1.f / theta);
                    // } else {

                    rhs[idx] += u_weights(i, j) * grid.velXHalfIndexed(i, j) / grid.deltaX()
                                + (1 - u_weights(i, j)) * u_solid / grid.deltaX();
                    // }
                }
                // right neighbour
                if (i + 1 < grid.sizeX()) {
                    // if (grid.cellType(i + 1, j) != SOLID) {
                    if (phi_i_plus_1_j < 0) {
                        alpha_matrix.insert(idx, fluid_indices(i_plus_1, j)) = -u_weights(i_plus_1, j) * scale_x;
                    }
                    theta = max(0.01f, fraction_inside(phi_i_j, phi_i_plus_1_j));
                    diagonal += u_weights(i_plus_1, j) * scale_x * (1.f / theta);
                    // } else {
                    rhs[idx] -= u_weights(i_plus_1, j) * grid.velXHalfIndexed(i_plus_1, j) / grid.deltaX()
                                + (1 - u_weights(i_plus_1, j)) * u_solid / grid.deltaX();
                    // }
                }
                // bottom neighbour
                // if (grid.cellType(i, j - 1) != SOLID) {
                if (j - 1 >= 0) {
                    if (phi_i_j_minus_1 < 0) {
                        alpha_matrix.insert(idx, fluid_indices(i, j_minus_1)) = -v_weights(i, j) * scale_y;
                    }
                    theta = max(0.01f, fraction_inside(phi_i_j_minus_1, phi_i_j));
                    diagonal += v_weights(i, j) * scale_y * (1.f / theta);
                    // } else {
                    rhs[idx] += v_weights(i, j) * grid.velYHalfIndexed(i, j) / grid.deltaY()
                                + (1 - v_weights(i, j)) * v_solid / grid.deltaY();
                    // }
                }
                // top neighbour
                // if (grid.cellType(i, j + 1) != SOLID) {
                if (j + 1 < grid.sizeY()) {
                    if (phi_i_j_plus_1 < 0) {
                        alpha_matrix.insert(idx, fluid_indices(i, j_plus_1)) = -v_weights(i, j_plus_1) * scale_y;
                    }
                    theta = max(0.01f, fraction_inside(phi_i_j, phi_i_j_plus_1));
                    diagonal += v_weights(i, j_plus_1) * scale_y * (1.f / theta);
                    // } else {
                    rhs[idx] -= v_weights(i, j_plus_1) * grid.velYHalfIndexed(i, j_plus_1) / grid.deltaY()
                                + (1 - v_weights(i, j_plus_1)) * v_solid / grid.deltaY();
                }
                alpha_matrix.insert(idx, idx) = diagonal;
            }
            // }
        }
    }
    // Solid boundary fix. THIS ASSUMES THAT ALL LIQUID CELLS ARE CONNECTED WHICH OF COURSE MIGHT BE WRONG FOR
    // VARIOUS TOPOLOGIES
    double mean = 0;
    int counter = 0;

    int fluid_cells = 0;
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (domain.fluidLevelSet()(i, j) < 0) fluid_cells++;
        }
    }
    if (fluid_cells == (grid.sizeX() - 2) * (grid.sizeY() - 2)) {
        for (int j = 0; j < grid.sizeY(); ++j) {
            for (int i = 0; i < grid.sizeX(); ++i) {
                int idx = fluid_indices(i, j); // i + j * grid.sizeX();
                int i_minus_1 = clamp(i - 1, 0, grid.sizeX() - 1);
                int j_minus_1 = clamp(j - 1, 0, grid.sizeY() - 1);

                int i_plus_1 = clamp(i + 1, 0, grid.sizeX() - 1);
                int j_plus_1 = clamp(j + 1, 0, grid.sizeY() - 1);

                if (domain.fluidLevelSet()(i, j) < 0) {
                    // Check if not AIR cell

                    if ((domain.fluidLevelSet()(i_minus_1, j) < 0 || domain.solidLevelSet()(i_minus_1, j) <= 0)
                        && (domain.fluidLevelSet()(i_plus_1, j) < 0 || domain.solidLevelSet()(i_plus_1, j) <= 0)
                        && (domain.fluidLevelSet()(i, j_minus_1) < 0 || domain.solidLevelSet()(i, j_minus_1) <= 0)
                        && (domain.fluidLevelSet()(i, i_plus_1) < 0 || domain.solidLevelSet()(i, i_plus_1) <= 0)) {
                        mean += rhs[idx];
                        counter++;
                    }
                }
            }
        }
        mean = mean / counter;
        for (int j = 0; j < grid.sizeY(); ++j) {
            for (int i = 0; i < grid.sizeX(); ++i) {
                int idx = fluid_indices(i, j); // i + j * grid.sizeX();
                int i_minus_1 = clamp(i - 1, 0, grid.sizeX() - 1);
                int j_minus_1 = clamp(j - 1, 0, grid.sizeY() - 1);

                int i_plus_1 = clamp(i + 1, 0, grid.sizeX() - 1);
                int j_plus_1 = clamp(j + 1, 0, grid.sizeY() - 1);
                if (domain.fluidLevelSet()(i, j) < 0) {
                    if ((domain.fluidLevelSet()(i_minus_1, j) < 0 || domain.solidLevelSet()(i_minus_1, j) <= 0)
                        && (domain.fluidLevelSet()(i_plus_1, j) < 0 || domain.solidLevelSet()(i_plus_1, j) <= 0)
                        && (domain.fluidLevelSet()(i, j_minus_1) < 0 || domain.solidLevelSet()(i, j_minus_1) <= 0)
                        && (domain.fluidLevelSet()(i, j_plus_1) < 0 || domain.solidLevelSet()(i, j_plus_1) <= 0)) {
                        rhs[idx] -= mean;
                    }
                }
            }
        }
    }

    Eigen::VectorXf pressure_grid(system_size);
    solver.compute(alpha_matrix);
    pressure_grid = solver.solve(rhs);

    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            int i_minus1 = clamp(i - 1, 0, grid.sizeX() - 1);
            int j_minus1 = clamp(j - 1, 0, grid.sizeY() - 1);
            int idx = fluid_indices(i, j);
            int idx_i_minus1 = fluid_indices(i_minus1, j);
            int idx_j_minus1 = fluid_indices(i, j_minus1);
            float p = idx >= 0 ? pressure_grid[idx] : 0;
            float p_i_minus1 = idx_i_minus1 >= 0 ? pressure_grid[idx_i_minus1] : 0;
            float p_j_minus1 = idx_j_minus1 >= 0 ? pressure_grid[idx_j_minus1] : 0;

            int i_minus_1 = clamp(i - 1, 0, grid.sizeX() - 1);
            int j_minus_1 = clamp(j - 1, 0, grid.sizeY() - 1);

            // if (idx >= 0 || idx_i_minus1 >= 0 || idx_j_minus1 >= 0) {
            if (domain.fluidLevelSet()(i, j) < 0 || domain.fluidLevelSet()(i_minus_1, j) < 0) {
                if (domain.solidLevelSet()(i, j) <= 0 || domain.solidLevelSet()(i_minus_1, j) <= 0) {
                    grid.setVelXBackBufferHalfIndexed(i, j, 0);
                } else {
                    float new_vel_x =
                        grid.velXHalfIndexed(i, j) - dt / FluidDomain::density * (p - p_i_minus1) / grid.deltaX();
                    grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);
                }
            }
            if (domain.fluidLevelSet()(i, j) < 0 || domain.fluidLevelSet()(i, j_minus_1) < 0) {
                if (domain.solidLevelSet()(i, j) <= 0 || domain.solidLevelSet()(i, j_minus_1) <= 0) {
                    grid.setVelYBackBufferHalfIndexed(i, j, 0);
                } else {
                    float new_vel_y =
                        grid.velYHalfIndexed(i, j) - dt / FluidDomain::density * (p - p_j_minus1) / grid.deltaY();
                    grid.setVelYBackBufferHalfIndexed(i, j, new_vel_y);
                }
            }
        }
    }
    grid.swapVelocityBuffers();
}

void FluidSimulator::enforceDirichlet(FluidDomain &domain) {
    auto &grid = domain.grid();
    // X vel
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            int i_minus1 = clamp(i - 1, 0, grid.sizeX() - 1);
            int j_minus1 = clamp(j - 1, 0, grid.sizeY() - 1);
            // X velocity
            if ((domain.solidLevelSet()(i_minus1, j) <= 0 && grid.velXHalfIndexed(i, j) < 0)
                || (domain.solidLevelSet()(i, j) <= 0 && grid.velXHalfIndexed(i, j) > 0))
                grid.setVelXHalfIndexed(i, j, 0);

            // Y velocity
            if ((domain.solidLevelSet()(i, j_minus1) <= 0 && grid.velYHalfIndexed(i, j) < 0)
                || (domain.solidLevelSet()(i, j) <= 0 && grid.velYHalfIndexed(i, j) > 0))
                grid.setVelYHalfIndexed(i, j, 0);
        }
    }
}

float FluidSimulator::calculate_kernel_function(MacGrid &grid, float x, float y) {
    float val1 = calculate_quad_bspline(x / grid.deltaX());
    float val2 = calculate_quad_bspline(y / grid.deltaY());
    return val1 * val2;
}