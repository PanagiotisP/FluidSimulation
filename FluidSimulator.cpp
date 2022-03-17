#pragma once
#include "FluidSimulator.h"
#include "vec.h"
#include "array2_utils.h"
//#define PRINT
#define VERTICAL_PLANE false
const float epsilon = 10e-37;

float liquid_bound_x0, liquid_bound_x1;
float liquid_bound_y0, liquid_bound_y1;

float solid_bound_x0, solid_bound_x1;
float solid_bound_y0, solid_bound_y1;


float distance_from_axis_aligned_rect_box(const float x0, const float x1, const float y0, const float y1, const Vec2f& point) {
    //inside
    if (x0 < point[0] && point[0] < x1
        && y0 < point[1] && point[1] < y1) {
        return max(x0 - point[0], point[0] - x1, y0 - point[1], point[1] - y1);
    }
    //outside
    else {
        //closest point to bounding box (p,q)
        float p, q;
        if (point[0] < x0) p = x0;
        else if (x1 < point[0]) p = x1;
        else p = point[0];

        if (point[1] < y0) q = y0;
        else if (y1 < point[1]) q = y1;
        else q = point[1];
        return sqrt((sqr(point[0] - p)) + (sqr(point[1] - q)));
    }
}

FluidSimulator::FluidSimulator() {}

FluidSimulator::~FluidSimulator() {}

void FluidSimulator::classifyCells(MacGrid& grid, LevelSet& levelSet) {
    //// First reset types (set all cells to AIR)
    //grid.clearCellTypeBuffer();
    //Fixed topology for now
    grid.clearCellTypeBuffer();
    //for (int j = 1; j < grid.sizeY() - 1; ++j) {
    //    for (int i = 1; i < grid.sizeX() - 1; ++i) {
    //        grid.setCellType(i, j, LIQUID);
    //    }
    //}
}

void FluidSimulator::print_temperature_field(MacGrid& grid, const char* variable_name) {
#ifdef PRINT
    std::cout << variable_name;
    for (int j = grid.sizeY() - 1; j >= 0; --j) {
        std::cout << std::endl;
        for (int i = 0; i < grid.sizeX(); ++i) {
            std::cout << grid.temperature(i, j) - ZERO_CELCIUS << " ";
        }
    }
    std::cout << std::endl;
#endif // PRINT
}

void FluidSimulator::print_concentration_field(MacGrid& grid, const char* variable_name) {
#ifdef PRINT
    std::cout << variable_name;
    for (int j = grid.sizeY() - 1; j >= 0; --j) {
        std::cout << std::endl;
        for (int i = 0; i < grid.sizeX(); ++i) {
            std::cout << grid.concentration(i, j) << " ";
        }
    }
    std::cout << std::endl;
#endif // PRINT
}

void FluidSimulator::print_velocity_field(MacGrid& grid, const char* variable_name) {
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

void FluidSimulator::diffuse_scalar(MacGrid& grid, float dt, float diffuse_rate) {
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (grid.cellType(i, j) == LIQUID) {
                int i_minus1 = clamp(i - 1, 0, grid.sizeX() - 1);
                int j_minus1 = clamp(j - 1, 0, grid.sizeY() - 1);

                int i_plus1 = clamp(i + 1, 0, grid.sizeX() - 1);
                int j_plus1 = clamp(j + 1, 0, grid.sizeY() - 1);


                float neighbouring_values = (grid.concentration(i_minus1, j) + grid.concentration(i, j_minus1) + grid.concentration(i_plus1, j) + grid.concentration(i, j_plus1)) / 4;
                float next_value = (grid.concentration(i, j) + diffuse_rate * neighbouring_values) / (1 + diffuse_rate);

                grid.setConcentrationBackBuffer(i, j, next_value);
            }
        }
    }
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            grid.setConcentration(i, j, grid.concentrationBackBuffer(i, j));
        }
    }

}

void FluidSimulator::add_temperature_source(MacGrid& grid, const ScalarSource& src, float dt) {
    float current_temperature = grid.temperatureInterpolated(src.x, src.y);
    float temperature_increase = (1 - exp(-rate_t * dt)) * (src.value - current_temperature);
    grid.addToTemperatureInterpolated(src.x, src.y, temperature_increase);
}

void FluidSimulator::add_concentration_source(MacGrid& grid, const ScalarSource& src, float dt) {
    float current_concentration = grid.concentrationInterpolated(src.x, src.y);
    float concentration_increase = dt * (src.value - current_concentration);
    grid.addToConcentrationInterpolated(src.x, src.y, concentration_increase);
}

void FluidSimulator::add_velocity_source(MacGrid& grid, const VectorSource& src, float dt) {
    grid.addToVelXInterpolated(src.x, src.y + 0.5 * grid.deltaY(), src.val_x);
    grid.addToVelYInterpolated(src.x + 0.5 * grid.deltaX(), src.y, src.val_y);
}

template <class T, class VectorT>
void FluidSimulator::apply_sources(MacGrid& grid, VectorT* sources, float dt, void (FluidSimulator::* function)(MacGrid&, const T&, float)) {
    for (auto it = sources->begin(); it != sources->end(); ++it) {
        (this->*function)(grid, *it, dt);
    }
}

void FluidSimulator::advance(MacGrid& grid, float t_frame) {
    std::vector<ScalarSource>* temp_scalar_vec = nullptr;
    std::vector<VectorSource>* temp_vector_vec = nullptr;
    advance(grid, t_frame, temp_scalar_vec, temp_scalar_vec, temp_vector_vec, false, false, false);
}
void FluidSimulator::advance(MacGrid& grid, float t_frame, std::vector<ScalarSource>* temperatureSrcs, std::vector<ScalarSource>* concentrationSrcs, std::vector<VectorSource>* velocitySrcs,
    bool isTemperatureSrcActive, bool isConcentrationSrcActive, bool isVelocitySrcActive) {
    float t = 0;
    while (t < t_frame) {
        float timestep = compute_cfl(grid);
        if (timestep + t >= t_frame) {
            timestep = t_frame - t;
        }
        printf("Taking substep of size %f (to %0.3f%% of the frame)\n", timestep, 100 * (t + timestep) / t_frame);
        LevelSet levelSet(grid.sizeX(), grid.sizeY(), 1.0f, 1.0f);
        classifyCells(grid, levelSet);


        //print_temperature_field(grid, "Before Advection");

        print_concentration_field(grid, "Concentration field");
        advect(grid, timestep);
        print_concentration_field(grid, "Concentration field");

        if (isTemperatureSrcActive) {
            apply_sources<ScalarSource>(grid, temperatureSrcs, timestep, &FluidSimulator::add_temperature_source);
            grid.swapTemperatureBuffers();
        }

        if (isConcentrationSrcActive) {
            apply_sources<ScalarSource>(grid, concentrationSrcs, timestep, &FluidSimulator::add_concentration_source);
            grid.swapConcentrationBuffers();
        }

        print_concentration_field(grid, "Concentration field");
        if (isVelocitySrcActive) {
            apply_sources<VectorSource>(grid, velocitySrcs, timestep, &FluidSimulator::add_velocity_source);
            grid.swapVelocityBuffers();
        }



        diffuse_scalar(grid, timestep, diffuse_rate_concentration);

        //apply forces (e.g. gravity)
        add_forces(grid, timestep);
        print_concentration_field(grid, "Concentration field");
        //print_velocity_field(grid, "After forces");
        //print_temperature_field(grid, "After forces");


        project(grid, timestep);
        //print_velocity_field(grid, "After projection");

        enforceDirichlet(grid);

        ////update temperature and concentration
        //float rate_t = 10;
        //float rate_s = 1;
        //float target_temperature = 273 + 37;
        ////temperature(grid.sizeX() / 2, grid.sizeY() / 2) += (1 - exp(-rate_t * timestep)) * (target_temperature - temperature(grid.sizeX() / 2, grid.sizeY() / 2));
        //temperature(0, 0) += (1 - exp(-rate_t * timestep)) * (target_temperature - temperature(0, 0));
        //temperature(0, 1) += (1 - exp(-rate_t * timestep)) * (target_temperature - temperature(0, 1));
        //concentration(grid.sizeX() / 2, grid.sizeY() / 2) += clamp(rate_s * timestep, 0.f, 1.f);

        t += timestep;
    }
}

void FluidSimulator::get_advected_position(MacGrid& grid, float x_initial, float y_initial, float dt, float* x_result, float* y_result) {
    float x_mid = x_initial - 0.5 * dt * grid.velXInterpolated(x_initial, y_initial);
    float y_mid = y_initial - 0.5 * dt * grid.velYInterpolated(x_initial, y_initial);

    *x_result = x_initial - dt * grid.velXInterpolated(x_mid, y_mid);
    *y_result = y_initial - dt * grid.velYInterpolated(x_mid, y_mid);
    return;
}

void produce_advected_position(MacGrid& grid, float x_initial, float y_initial, float dt, float* x_result, float* y_result) {
    float vel_x_k1 = grid.velXInterpolated(x_initial, y_initial);
    float vel_y_k1 = grid.velYInterpolated(x_initial, y_initial);

    float vel_x_k2 = grid.velXInterpolated(x_initial + 0.5 * dt * vel_x_k1, y_initial + 0.5 * dt * vel_y_k1);
    float vel_y_k2 = grid.velYInterpolated(x_initial + 0.5 * dt * vel_x_k1, y_initial + 0.5 * dt * vel_y_k1);

    float vel_x_k3 = grid.velXInterpolated(x_initial + 0.75 * dt * vel_x_k2, y_initial + 0.75 * dt * vel_y_k2);
    float vel_y_k3 = grid.velYInterpolated(x_initial + 0.75 * dt * vel_x_k2, y_initial + 0.75 * dt * vel_y_k2);

    *x_result = x_initial + 2 / 9 * dt * vel_x_k1 + 3 / 9 * dt * vel_x_k2 + 4 / 9 * dt * vel_x_k3;
    *y_result = y_initial + 2 / 9 * dt * vel_y_k1 + 3 / 9 * dt * vel_y_k2 + 4 / 9 * dt * vel_y_k3;
    return;
}

//void FluidSimulator::advectLevelSet(MacGrid& grid, LevelSet& levelSet, float dt) {
//    LevelSet newlevelSet(
//        levelSet.sizeX(),
//        levelSet.sizeY(),
//        levelSet.lengthX(),
//        levelSet.lengthY());
//    for (int j = 0; j < levelSet.sizeX(); ++j) 	{
//        for (int i = 0; i < levelSet.sizeY(); ++i) 		{
//            float vel_x = grid.velX(i, j);
//            float vel_y = grid.velY(i, j);
//            float grad_x = levelSet.computeUpwindGradientX(i, j, vel_x);
//            float grad_y = levelSet.computeUpwindGradientY(i, j, vel_y);
//
//            float change_rate = -(vel_x * grad_x + vel_y * grad_y);
//
//            // Forward Euler
//            newlevelSet(i, j) = levelSet(i, j) + change_rate * dt * 5;
//        }
//    }
//    for (int j = 0; j < levelSet.sizeX(); ++j) {
//        for (int i = 0; i < levelSet.sizeY(); ++i) {
//            levelSet(i, j) = newlevelSet(i, j);
//        }
//    }
//}

float FluidSimulator::compute_cfl(MacGrid& grid) {
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
void FluidSimulator::advect_particles(MacGrid& grid, float dt) {
    //float substep = 0;
    //float t = 0;
    //while (substep + t <= dt) {

    //    if (2 * substep + t >= dt) {
    //        substep = 0.5 * (dt - t);
    //    }
    //}


}
void FluidSimulator::advect(MacGrid& grid, float dt) {
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            grid.setVelXBackBufferHalfIndexed(i, j, grid.velXHalfIndexed(i, j));
            grid.setVelYBackBufferHalfIndexed(i, j, grid.velYHalfIndexed(i, j));
            grid.setTemperatureBackBuffer(i, j, ambient_temp);
            grid.setConcentrationBackBuffer(i, j, 0);
        }
    }

    // Hypothetical particle for x-component
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (grid.cellType(i, j) == LIQUID || grid.cellType(i - 1, j) == LIQUID) {
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
            if (grid.cellType(i, j) == LIQUID || grid.cellType(i - 1, j) == LIQUID) {
                float x_current = (i + 0.5) * grid.deltaX();
                float y_current = j * grid.deltaY();
                float x_previous, y_previous;

                get_advected_position(grid, x_current, y_current, dt, &x_previous, &y_previous);

                float v_y = grid.velYInterpolated(x_previous, y_previous);
                grid.setVelYBackBufferHalfIndexed(i, j, v_y);
            }
        }
    }

    //Advect other quantities like Temperature and Concentration
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            if (grid.cellType(i, j) == LIQUID || grid.cellType(i - 1, j) == LIQUID) {
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

    //Write results to front buffer
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            grid.setVelXHalfIndexed(i, j, grid.velXBackBufferHalfIndexed(i, j));
            grid.setVelYHalfIndexed(i, j, grid.velYBackBufferHalfIndexed(i, j));

            //Advect scalar field too
            grid.setTemperature(i, j, grid.temperatureBackBuffer(i, j));
            grid.setConcentration(i, j, grid.concentrationBackBuffer(i, j));
        }
    }
}

void FluidSimulator::add_forces(MacGrid& grid, float dt) {
    //gravity
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            // Only add force to the liquid cells
            if (grid.cellType(i, j) == LIQUID) {
                if (VERTICAL_PLANE) {
                    grid.setVelYHalfIndexed(i, j,
                        grid.velYHalfIndexed(i, j) - grav * dt);
                }
            }
        }
    }
}

void FluidSimulator::project(MacGrid& grid, float dt) {
    const int system_size = grid.sizeX() * grid.sizeY();

    //Pressure solver
    std::vector<double> rhs(system_size);
    std::vector<double> pressure_grid(system_size);
    SparseMatrixd alpha_matrix(system_size);
    PCGSolver<double> solver;

    alpha_matrix.zero();
    rhs.assign(rhs.size(), 0);
    pressure_grid.assign(pressure_grid.size(), 0);


    float scale_x = dt / (density * sqr(grid.deltaX()));
    float scale_y = dt / (density * sqr(grid.deltaY()));
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            int idx = i + j * grid.sizeX();

            rhs[idx] = 0;
            pressure_grid[idx] = 0;
            //solid velocity zero for now
            if (grid.cellType(i, j) == LIQUID) {
                //left neighbour
                if (grid.cellType(i - 1, j) != SOLID) {
                    if (grid.cellType(i - 1, j) == LIQUID) {
                        alpha_matrix.add_to_element(idx, idx - 1, -scale_x);
                    }
                    alpha_matrix.add_to_element(idx, idx, scale_x);
                }
                else {
                    rhs[idx] -= grid.velXHalfIndexed(i, j) / grid.deltaX();
                }

                //right neighbour
                if (grid.cellType(i + 1, j) != SOLID) {
                    if (grid.cellType(i + 1, j) == LIQUID) {
                        alpha_matrix.add_to_element(idx, idx + 1, -scale_x);
                    }
                    alpha_matrix.add_to_element(idx, idx, scale_x);
                }
                else {
                    rhs[idx] += grid.velXHalfIndexed(i + 1, j) / grid.deltaX();
                }

                //bottom neighbour
                if (grid.cellType(i, j - 1) != SOLID) {
                    if (grid.cellType(i, j - 1) == LIQUID) {
                        alpha_matrix.add_to_element(idx, idx - grid.sizeX(), -scale_y);
                    }
                    alpha_matrix.add_to_element(idx, idx, scale_y);
                }
                else {
                    rhs[idx] -= grid.velYHalfIndexed(i, j) / grid.deltaY();
                }

                //top neighbour
                if (grid.cellType(i, j + 1) != SOLID) {
                    if (grid.cellType(i, j + 1) == LIQUID) {
                        alpha_matrix.add_to_element(idx, idx + grid.sizeX(), -scale_y);
                    }
                    alpha_matrix.add_to_element(idx, idx, scale_y);
                }
                else {
                    rhs[idx] += grid.velYHalfIndexed(i, j + 1) / grid.deltaY();
                }

                // Set right hand side equal to negative divirgence
                rhs[idx] -= (grid.divVelX(i, j) + grid.divVelY(i, j));
            }
        }
    }
    //Solve the system using Robert Bridson's incomplete Cholesky PCG solver

    double tolerance;
    int iterations;

    // Solid boundary fix. THIS ASSUMES THAT ALL LIQUID CELLS ARE CONNECTED WHICH OF COURSE MIGHT BE WRONG FOR VARIOUS TOPOLOGIES
    double mean = 0;
    int counter = 0;
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            int idx = i + j * grid.sizeX();
            if (grid.cellType(i, j) == LIQUID) {
                if (grid.cellType(i - 1, j) != AIR && grid.cellType(i + 1, j) != AIR && grid.cellType(i, j - 1) != AIR && grid.cellType(i - 1, j + 1) != AIR) {
                    mean += rhs[idx];
                    counter++;
                }
            }
        }
    }
    mean = mean / counter;
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            int idx = i + j * grid.sizeX();
            if (grid.cellType(i, j) == LIQUID) {
                if (grid.cellType(i - 1, j) != AIR && grid.cellType(i + 1, j) != AIR && grid.cellType(i, j - 1) != AIR && grid.cellType(i - 1, j + 1) != AIR) {
                    rhs[idx] -= mean;
                }
            }
        }
    }

    solver.set_solver_parameters(1e-15, 1000);
    bool success = solver.solve(alpha_matrix, rhs, pressure_grid, tolerance, iterations);
#ifdef PRINT
    printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
#endif // PRINT
    if (!success) {
        printf("WARNING: Pressure solve failed!************************************************\n");
    }

    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            //int i_minus1 = clamp(i - 1, 0, grid.sizeX() - 1);
            //int j_minus1 = clamp(j - 1, 0, grid.sizeY() - 1);
            int idx = i + j * grid.sizeX();
            int idx_i_minus1 = idx - 1;
            int idx_j_minus1 = idx - grid.sizeX();
            float p = idx >= 0 ? pressure_grid[idx] : 0;
            float p_i_minus1 = idx_i_minus1 >= 0 ? pressure_grid[idx_i_minus1] : 0;
            float p_j_minus1 = idx_j_minus1 >= 0 ? pressure_grid[idx_j_minus1] : 0;


            //if (idx >= 0 || idx_i_minus1 >= 0 || idx_j_minus1 >= 0) {
            if (grid.cellType(i, j) == LIQUID || grid.cellType(i - 1, j) == LIQUID) {
                if (grid.cellType(i, j) == SOLID || grid.cellType(i - 1, j) == SOLID) {
                    grid.setVelXBackBufferHalfIndexed(i, j, 0);
                }
                else {
                    float new_vel_x = grid.velXHalfIndexed(i, j) - dt / density * (p - p_i_minus1) / grid.deltaX();
                    grid.setVelXBackBufferHalfIndexed(i, j, new_vel_x);
                }
            }
            if (grid.cellType(i, j) == LIQUID || grid.cellType(i, j - 1) == LIQUID) {
                if (grid.cellType(i, j) == SOLID || grid.cellType(i, j - 1) == SOLID) {
                    grid.setVelYBackBufferHalfIndexed(i, j, 0);
                }
                else {
                    float new_vel_y = grid.velYHalfIndexed(i, j) - dt / density * (p - p_j_minus1) / grid.deltaY();
                    grid.setVelYBackBufferHalfIndexed(i, j, new_vel_y);
                }
            }
        }
    }
    grid.swapVelocityBuffers();
}

void FluidSimulator::enforceDirichlet(MacGrid& grid) {
    // X vel
    for (int j = 0; j < grid.sizeY(); ++j) {
        for (int i = 0; i < grid.sizeX(); ++i) {
            int i_minus1 = clamp(i - 1, 0, grid.sizeX() - 1);
            int j_minus1 = clamp(j - 1, 0, grid.sizeY() - 1);
            // X velocity
            if ((grid.cellType(i_minus1, j) == SOLID &&
                grid.velXHalfIndexed(i, j) < 0) ||
                (grid.cellType(i, j) == SOLID &&
                    grid.velXHalfIndexed(i, j) > 0))
                grid.setVelXHalfIndexed(i, j, 0);

            // Y velocity
            if ((grid.cellType(i, j_minus1) == SOLID &&
                grid.velYHalfIndexed(i, j) < 0) ||
                (grid.cellType(i, j) == SOLID &&
                    grid.velYHalfIndexed(i, j) > 0))
                grid.setVelYHalfIndexed(i, j, 0);
        }
    }
}

