#pragma once
#include "FluidSimulator.h"

#include "FluidDomain.h"
#include "Functors.h"
#include "util.h"
#include "vec.h"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCore>
#include <chrono>
#include <openvdb/math/Stats.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/Count.h>
#include <openvdb/tools/FastSweeping.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/Morphology.h>
#include <openvdb/tools/ParticleAtlas.h>
#include <openvdb/tools/PointsToMask.h>
#include <openvdb/tools/TopologyToLevelSet.h>
#include <random>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

// #define PRINT
#define VERTICAL_PLANE true
#define GATHER_TRANSFER
const float epsilon = 10e-37;

float liquid_bound_x0, liquid_bound_x1;
float liquid_bound_y0, liquid_bound_y1;

float solid_bound_x0, solid_bound_x1;
float solid_bound_y0, solid_bound_y1;

using namespace std::chrono;

FluidSimulator::FluidSimulator(openvdb::CoordBBox bbox): _bbox(bbox) {}

FluidSimulator::~FluidSimulator() {}

void FluidSimulator::print_velocity_field(MacGrid &grid, const char *variable_name) {
#ifdef PRINT
    auto accessor = grid.velFront()->getAccessor();
    auto bbox = grid.velFront()->evalActiveVoxelBoundingBox();
    std::cout << variable_name;
    // assert(bbox.min() < bbox.max());
    // for (int i = _bbox.min()[0]; i < _bbox.max()[0]; ++i) {
    for (int i = 14; i < 17; ++i) {
        std::cout << "\n" << i << std::endl;
        // for (int j = _bbox.max()[1]; j >= _bbox.min()[1]; --j) {
        for (int j = 17; j >= 13; --j) {
            std::cout << std::endl;
            // for (int k = _bbox.min()[2]; k < _bbox.max()[2]; ++k) {
            for (int k = 15; k < 18; ++k) {
                std::cout << "[" << grid.velHalfIndexed(accessor, i, j, k)[0] << ", "
                          << grid.velHalfIndexed(accessor, i, j, k)[1] << ", "
                          << grid.velHalfIndexed(accessor, i, j, k)[2] << "]"
                          << " ";
            }
        }
    }
    std::cout << std::endl;
#endif // PRINT
}

void FluidSimulator::transfer_from_particles_to_grid(FluidDomain &domain) {
    MacGrid &grid = domain.grid();

    grid.velBack()->clear();

#ifdef GATHER_TRANSFER
    openvdb::Vec3dGrid::Ptr vel_weight(openvdb::Vec3DGrid::create(openvdb::Vec3d(epsilon)));
    vel_weight->setTransform(domain._i2w_transform);

    // Zero out radius for the atlas creation as we care only for the particles whose center lies inside a voxel
    domain.particleSet().setZeroRadius(true);

    auto p_atlas = openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::create<ParticleSet>(
        domain.particleSet(), domain.voxelSize());

    // Create a point mask to pre-allocate back velocity and weight grids, in order to perform parallel writes.
    auto pointMask = openvdb::tools::createPointMask(domain.particleSet(), grid.velBack()->transform());
    // Dilate the point mask to cover the particles' contribution, as determined by the kernel function
    openvdb::tools::dilateActiveValues(pointMask->tree(), 2, openvdb::tools::NN_FACE_EDGE_VERTEX);

    grid.velBack()->tree().topologyUnion(pointMask->tree());
    grid.velBack()->tree().voxelizeActiveTiles();

    vel_weight->tree().topologyUnion(pointMask->tree());
    vel_weight->tree().voxelizeActiveTiles();

    if (pointMask->activeVoxelCount() != 0) {
        auto bbox = pointMask->evalActiveVoxelBoundingBox();
        tbb::parallel_for(tbb::blocked_range<int>(bbox.min()[0], bbox.max()[0], 1),
                          GatherTransfer(domain, p_atlas, bbox, grid.velBack(), vel_weight));
    }

    // Re-enable particles' radii
    domain.particleSet().setZeroRadius(false);

    // Normalize velocity grid with the weight values
    openvdb::tools::compDiv<openvdb::Vec3dGrid>(*grid.velBack(), *vel_weight);

#else
    ShootingTransfer f(domain, domain.particleSet());
    tbb::parallel_reduce(tbb::blocked_range<int>(0, domain.particleSet().size(), 1), f);
    openvdb::tools::compDiv<openvdb::Vec3dGrid>(*f.vel, *f.vel_weight);
    grid.setVelBack(f.vel->deepCopy());
#endif

    grid.swapVelocityBuffers();
}

void FluidSimulator::transfer_from_grid_to_particles(FluidDomain &domain, float flip_pic_ratio = 0.98) {
    MacGrid &grid = domain.grid();

    tbb::parallel_for(tbb::blocked_range<int>(0, domain.particleSet().size(), 1), [&](tbb::blocked_range<int> range) {
        MacGrid::StaggeredSampler vel_sampler(grid.velFront()->getAccessor(), grid.velFront()->transform());
        MacGrid::StaggeredSampler vel_diff_sampler(grid.velDiff()->getAccessor(), grid.velDiff()->transform());

        MacGrid::Sampler displacement_sampler(grid.displacementGrid()->getAccessor(),
                                              grid.displacementGrid()->transform());

        for (int i = range.begin(); i < range.end(); ++i) {
            auto &p = domain.particleSet()[i];
            openvdb::Vec3d pic_vel = vel_sampler.isSample(p.pos());
            openvdb::Vec3d flip_vel = p.vel() + vel_diff_sampler.isSample(p.pos());

            p.setVelocity(pic_vel * (1 - flip_pic_ratio) + flip_vel * flip_pic_ratio);
        }
    });
}

void FluidSimulator::reseeding(FluidDomain &domain) {
    domain.particleSet().setZeroRadius(true);
    MacGrid &grid(domain.grid());

    openvdb::MaskGrid::Ptr active_mask =
        openvdb::tools::createPointMask(domain.particleSet(), domain.fluidLevelSet().getLevelSet()->transform());

    // Perform opening
    // openvdb::tools::erodeActiveValues(active_mask->tree(), 1);
    // openvdb::tools::dilateActiveValues(active_mask->tree(), 1);

    auto p_atlas = openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::create<ParticleSet>(
        domain.particleSet(), domain.voxelSize());

    ReseedingFunctor f(domain, p_atlas, active_mask, grid.velFront(), 4, 8);
    openvdb::tree::IteratorRange<openvdb::MaskGrid::ValueOnCIter> range(active_mask->cbeginValueOn());
    tbb::parallel_reduce(range, f);

    std::sort(f.particles_to_be_deleted.begin(), f.particles_to_be_deleted.end());
    for (auto it = f.particles_to_be_deleted.rbegin(); it != f.particles_to_be_deleted.rend(); ++it) {
        domain.particleSet().removeParticle(domain.particleSet().begin() + *it);
    }
    for (auto it = f.particles_to_be_added.begin(); it != f.particles_to_be_added.end(); ++it) {
        domain.particleSet().addParticle(*it);
    }
    domain.particleSet().setZeroRadius(false);
}

void FluidSimulator::extrapolate_data(FluidDomain &domain, int iterations_n) {
    auto &grid = domain.grid();
    auto fluidAccessor = domain.fluidLevelSet().getAccessor();

    for (int iter = 0; iter < iterations_n; ++iter) {
        // Resetting accessors due to the swap happening at the end of this iteration
        grid.setVelBack(grid.velFront()->deepCopy());

        grid.setValidUBack(grid.validUFront()->deepCopy());
        grid.setValidVBack(grid.validVFront()->deepCopy());
        grid.setValidWBack(grid.validWFront()->deepCopy());

        auto vel_front_accessor = domain.grid().velFront()->getAccessor();
        auto vel_back_accessor = domain.grid().velBack()->getAccessor();

        auto valid_mask_u_front_buffer_accessor = grid.validUFront()->getAccessor();
        auto valid_mask_u_back_buffer_accessor = grid.validUBack()->getAccessor();
        auto valid_mask_v_front_buffer_accessor = grid.validVFront()->getAccessor();
        auto valid_mask_v_back_buffer_accessor = grid.validVBack()->getAccessor();
        auto valid_mask_w_front_buffer_accessor = grid.validWFront()->getAccessor();
        auto valid_mask_w_back_buffer_accessor = grid.validWBack()->getAccessor();

        auto solidAccessor = domain.solidLevelSet().getAccessor();

        auto bbox = domain.solidLevelSet().getLevelSet()->evalActiveVoxelBoundingBox();

        // Lambda function that performs an update for a single cell, for a single dimension
        auto updateOneDimension = [&](openvdb::MaskGrid::Accessor valid_mask_front_accessor,
                                      openvdb::MaskGrid::Accessor valid_mask_back_accessor, openvdb::Coord coord,
                                      int idx) {
            if (!valid_mask_front_accessor.isValueOn(coord)) {
                float new_vel = 0;
                int n_valid_neighbors = 0;

                // Get values of all valid neighbors
                if (valid_mask_front_accessor.isValueOn(coord.offsetBy(-1, 0, 0))) {
                    new_vel += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(-1, 0, 0))[idx];
                    n_valid_neighbors++;
                }
                if (valid_mask_front_accessor.isValueOn(coord.offsetBy(0, -1, 0))) {
                    new_vel += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, -1, 0))[idx];
                    n_valid_neighbors++;
                }
                if (valid_mask_front_accessor.isValueOn(coord.offsetBy(0, 0, -1))) {
                    new_vel += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 0, -1))[idx];
                    n_valid_neighbors++;
                }
                if (valid_mask_front_accessor.isValueOn(coord.offsetBy(1, 0, 0))) {
                    new_vel += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(1, 0, 0))[idx];
                    n_valid_neighbors++;
                }
                if (valid_mask_front_accessor.isValueOn(coord.offsetBy(0, 1, 0))) {
                    new_vel += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 1, 0))[idx];
                    n_valid_neighbors++;
                }
                if (valid_mask_front_accessor.isValueOn(coord.offsetBy(0, 0, 1))) {
                    new_vel += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 0, 1))[idx];
                    n_valid_neighbors++;
                }

                // Average the value for the current cell
                if (n_valid_neighbors > 0) {
                    grid.setVelOneDimHalfIndexed(vel_back_accessor, coord, new_vel / n_valid_neighbors, idx);
                    valid_mask_back_accessor.setValueOn(coord);
                }
            }
        };


        for (auto it = bbox.beginZYX(); it != bbox.endZYX(); ++it) {
            auto coord = (*it);

            // Extrapolate on x dimension
            updateOneDimension(valid_mask_u_front_buffer_accessor, valid_mask_u_back_buffer_accessor, coord, 0);
            // Extrapolate on y dimension
            updateOneDimension(valid_mask_v_front_buffer_accessor, valid_mask_v_back_buffer_accessor, coord, 1);
            // Extrapolate on z dimension
            updateOneDimension(valid_mask_w_front_buffer_accessor, valid_mask_w_back_buffer_accessor, coord, 2);
        }

        // Swap valid mask
        grid.swapValidMasks();
        grid.swapVelocityBuffers();
    }
}

void FluidSimulator::index_fluid_cells(FluidDomain &domain, openvdb::Int32Grid::Ptr fluid_indices) {
    auto &grid = domain.grid();
    // Basically a 3d to linear converter for fluid cells only
    auto fluid_indices_accessor = fluid_indices->getAccessor();

    auto fluidAccessor = domain.fluidLevelSet().getAccessor();
    auto solidAccessor = domain.solidLevelSet().getAccessor();
    auto vel_front_accessor = grid.velFront()->getAccessor();
    auto vel_back_accessor = grid.velBack()->getAccessor();

    auto u_weights_accessor = grid.uWeights()->getAccessor();
    auto v_weights_accessor = grid.vWeights()->getAccessor();
    auto w_weights_accessor = grid.wWeights()->getAccessor();

    auto u_fluid_weights_accessor = grid.uFluidWeights()->getAccessor();
    auto v_fluid_weights_accessor = grid.vFluidWeights()->getAccessor();
    auto w_fluid_weights_accessor = grid.wFluidWeights()->getAccessor();

    auto border_mask = openvdb::tools::createPointMask(domain.particleSet(), *domain._i2w_transform);
    openvdb::tools::dilateActiveValues(border_mask->tree(), 1, openvdb::tools::NN_FACE_EDGE_VERTEX,
                                       openvdb::tools::TilePolicy::EXPAND_TILES);

    openvdb::tree::IteratorRange<openvdb::MaskGrid::ValueOnCIter> fluidRange(border_mask->cbeginValueOn());
    int index = 0;
    if (domain.fluidLevelSet().getLevelSet()->activeVoxelCount() != 0) {
        for (; fluidRange; ++fluidRange) {
            auto coord = fluidRange.iterator().getCoord();
            // Index all cells with a computable pressure, that is cells with a non zero face area fraction
            // and that are either inside the fluid or a solid (for finite volume method)
            auto u_fluid_frac = (1 - u_weights_accessor.getValue(coord));
            auto v_fluid_frac = (1 - v_weights_accessor.getValue(coord));
            auto w_fluid_frac = (1 - w_weights_accessor.getValue(coord));


            if (fluidAccessor.getValue(coord) < 0
                || (((0 < u_fluid_frac && u_fluid_frac < 1) || (0 < v_fluid_frac && v_fluid_frac < 1)
                     || (0 < w_fluid_frac && w_fluid_frac < 1))
                    && solidAccessor.getValue(coord) <= 0)) {
                fluid_indices_accessor.setValue(coord, index++);
            }
        }
    }
}

void FluidSimulator::advance_flip_pic(FluidDomain &domain, float t_frame, float flip_pic_ratio = 0.98) {
    float t = 0;
    while (t < t_frame) {
        float timestep = compute_cfl(domain);
        if (timestep + t >= t_frame) {
            timestep = t_frame - t;
        }
        printf("Taking substep of size %f (to %0.3f%% of the frame)\n", timestep, 100 * (t + timestep) / t_frame);

        // Fluid sources update
        auto t_start = high_resolution_clock::now();
        domain.update(timestep);
        auto t_end = high_resolution_clock::now();
        auto update_duration = duration_cast<milliseconds>(t_end - t_start);

        // Fluid level set construction
        t_start = high_resolution_clock::now();
        domain.constructFluidLevelSetFromMarkerParticles();
        t_end = high_resolution_clock::now();
        auto constructFluidLevelSetFromMarkerParticles_duration =
            duration_cast<milliseconds>(t_end - t_start);

        // Reseeding (should probably not be called every frame but with a set interval)
        t_start = high_resolution_clock::now();
        // reseeding(domain);
        t_end = high_resolution_clock::now();
        auto reseeding_duration = duration_cast<milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "before transfer");

        // Transfer information from particles to grid
        t_start = high_resolution_clock::now();
        transfer_from_particles_to_grid(domain);
        t_end = high_resolution_clock::now();
        auto transfer_from_particles_to_grid_duration =
            duration_cast<milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after transfer");

        // Pre-allocate valid grids for the extrapolation function
        domain.grid().validUFront()->tree().topologyUnion(domain.grid().velFront()->tree());
        domain.grid().validVFront()->tree().topologyUnion(domain.grid().velFront()->tree());
        domain.grid().validWFront()->tree().topologyUnion(domain.grid().velFront()->tree());

        // Extrapolate and costrain velocities
        t_start = high_resolution_clock::now();
        extrapolate_data(domain, 2);
        t_end = high_resolution_clock::now();
        auto extrapolate_data_duration = duration_cast<milliseconds>(t_end - t_start);

        print_velocity_field(domain.grid(), "after extrapolation");

        t_start = high_resolution_clock::now();
        // Constrain on both solid boundary and solid objects
        constrain_velocity(domain);
        for (auto solidObj : domain.solidObjs()) {
            solidObj->generateVelocityField();
            constrain_velocity(domain, *solidObj);
        }
        t_end = high_resolution_clock::now();
        auto constrain_velocity_duration = duration_cast<milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after solid constrain");

        // Store previous velocities for FLIP update (which uses the difference to update velocities)
        domain.grid().updatePreviousVelocityBuffer();

        // Eulerian grid part START
        
        // Add body forices
        t_start = high_resolution_clock::now();
        add_forces(domain, timestep);
        t_end = high_resolution_clock::now();
        auto add_forces_duration = duration_cast<milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after forces");

        // Prepare grids necessary for solving the pressure poisson equation
        openvdb::Int32Grid::Ptr fluid_indices = openvdb::createGrid<openvdb::Int32Grid>(-1);
        prepare_pressure_solve(domain, fluid_indices);

        // Calculate pressure and project fluid velocities
        t_start = high_resolution_clock::now();
        solve_pressure_divirgence(domain, fluid_indices, timestep);
        t_end = high_resolution_clock::now();
        auto project_duration = duration_cast<milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after project");
        
        // Extrapolate and constrain velocities again
        t_start = high_resolution_clock::now();
        extrapolate_data(domain, 2);
        t_end = high_resolution_clock::now();
        extrapolate_data_duration += duration_cast<milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after extrapolation");

        t_start = high_resolution_clock::now();
        constrain_velocity(domain);
        for (auto solidObj : domain.solidObjs()) {
            solidObj->generateVelocityField();
            constrain_velocity(domain, *solidObj);
        }
        t_end = high_resolution_clock::now();
        constrain_velocity_duration += duration_cast<milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after solid constrain");

        // Eulerian grid part END

        // Update difference buffers for FLIP update
        domain.grid().updateDiffBuffers();

        // Perform FLIP/PIC update
        t_start = high_resolution_clock::now();
        transfer_from_grid_to_particles(domain, flip_pic_ratio);
        t_end = high_resolution_clock::now();
        auto transfer_from_grid_to_particles_duration =
            duration_cast<milliseconds>(t_end - t_start);

        // Advect particles and solids based on the calculated velocities
        t_start = high_resolution_clock::now();
        domain.advectParticles(timestep);
        for (auto solidObj : domain.solidObjs()) {
            solidObj->update_positions_orientations(timestep);
        }
        t_end = high_resolution_clock::now();
        auto advectParticles_duration = duration_cast<milliseconds>(t_end - t_start);

        printf(
            "update %.3f, construct %.3f, reseeding_duration %.3f, transfer grid %.3f, add forces %.3f, project %.3f,extrapolate %.3f, constrain_velocity_duration %.3f, transfer particles %.3f, advect particles %.3f \n",
            update_duration.count() / 1000.f, constructFluidLevelSetFromMarkerParticles_duration.count() / 1000.f,
            reseeding_duration.count() / 1000.f, transfer_from_particles_to_grid_duration.count() / 1000.f,
            add_forces_duration.count() / 1000.f, project_duration.count() / 1000.f,
            extrapolate_data_duration.count() / 1000.f, constrain_velocity_duration.count() / 1000.f,
            transfer_from_grid_to_particles_duration.count() / 1000.f, advectParticles_duration.count() / 1000.f);
        t += timestep;
    }
}

float FluidSimulator::compute_cfl(FluidDomain &domain) {
    auto &grid = domain.grid();
    
    openvdb::Vec3d minVector, maxVector;
    
    // Find max velocity
    grid.velFront()->evalMinMax(minVector, maxVector);
    return domain.voxelSize() / (max(maxVector[0], maxVector[1], maxVector[2]) + epsilon);
}

void FluidSimulator::add_forces(FluidDomain &domain, float dt) {
    auto &grid = domain.grid();
    auto accessor = domain.fluidLevelSet().getAccessor();
    auto vel_accessor = grid.velFront()->getAccessor();
    for (auto iter = grid.velFront()->cbeginValueOn(); iter; ++iter) {
        auto coord = iter.getCoord();
        // Only add force to the liquid cells
        // if (accessor.getValue(coord) < 0) {
        if (VERTICAL_PLANE) {
            grid.setVelOneDimHalfIndexed(vel_accessor, coord, grid.velHalfIndexed(vel_accessor, coord)[1] - grav * dt,
                                         1);
        }
        // }
    }

    // Add forces to solids
    // Only gravity for now
    for (auto solidObj : domain.solidObjs()) {
        solidObj->add_forces(dt);
    }
}

void FluidSimulator::prepare_pressure_solve(FluidDomain &domain, openvdb::Int32Grid::Ptr fluid_indices) {
    auto &grid = domain.grid();

    // Index all cells with a valid pressure
    index_fluid_cells(domain, fluid_indices);

    compute_face_fractions(domain);

}
void FluidSimulator::solve_pressure_divirgence(FluidDomain &domain, openvdb::Int32Grid::Ptr fluid_indices, float dt) {
    auto &grid = domain.grid();

    auto fluid_indices_accessor = fluid_indices->getAccessor();
    int system_size = fluid_indices->activeVoxelCount();

    if (system_size == 0) {
        return;
    }

    Eigen::VectorXd pressure_grid_divirgence_constraint(compute_pressure_divirgence_constraint(domain, fluid_indices_accessor, system_size, dt));

    project_divirgence_constraint(domain, fluid_indices_accessor, pressure_grid_divirgence_constraint, dt);
}
        }
    }

    if (system_size == 0) {
        return;
    }

    alpha_matrix.resize(system_size, system_size);
    alpha_matrix.reserve(Eigen::VectorXi::Constant(system_size, 7));
    
    for (auto solidObj : domain.solidObjs()) {
        solidObj->jMatrix().resize(6, system_size);
        solidObj->jMatrix().reserve(Eigen::VectorXi::Constant(6, system_size));
    }
    
    Eigen::VectorXd rhs(system_size);

void FluidSimulator::compute_face_fractions(FluidDomain &domain) {
    auto &grid = domain.grid();

    // Compute face area fractions for solid objects
    for (auto solidObj : domain.solidObjs()) {
        solidObj->compute_face_fractions();
    }

    grid.uWeights()->clear();
    grid.vWeights()->clear();
    grid.wWeights()->clear();

    grid.uFluidWeights()->clear();
    grid.vFluidWeights()->clear();
    grid.wFluidWeights()->clear();

    auto t_start = high_resolution_clock::now();

    // Compute face area fractions for solid boundary.
    // TODO Probably merge it with the one above

    // Resample solid level set to get values at bottom left corner of voxels
    // This is faster than sampling the bottom left for each cell separately

    auto computationLambda = [&](LevelSet &levelSet, openvdb::FloatGrid::Ptr uWeights, openvdb::FloatGrid::Ptr vWeights,
                                 openvdb::FloatGrid::Ptr wWeights) {
        // This resampling might be a very very bad low pass filter, destroying information
        auto resampledLevelSet = openvdb::createLevelSet<openvdb::FloatGrid>(domain.voxelSize());
        openvdb::tools::GridTransformer transformer(
            levelSet.getLevelSet()->transformPtr()->baseMap()->getAffineMap()->getMat4()
            * resampledLevelSet->transformPtr()->baseMap()->getAffineMap()->getMat4().inverse());

        transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(*levelSet.getLevelSet(),
                                                                                  *resampledLevelSet);
        // Preallocate weight grids for safe parallel writes
        auto border_mask = openvdb::tools::extractIsosurfaceMask(*domain.fluidLevelSet().getLevelSet(), 0);
        border_mask->tree().topologyUnion(
            openvdb::tools::extractEnclosedRegion(*domain.fluidLevelSet().getLevelSet(), 0)->tree());
        openvdb::tools::dilateActiveValues(border_mask->tree(), 1, openvdb::tools::NN_FACE_EDGE_VERTEX, openvdb::tools::TilePolicy::EXPAND_TILES);

        uWeights->tree().topologyUnion(border_mask->tree());
        vWeights->tree().topologyUnion(border_mask->tree());
        wWeights->tree().topologyUnion(border_mask->tree());

        uWeights->tree().voxelizeActiveTiles();
        vWeights->tree().voxelizeActiveTiles();
        wWeights->tree().voxelizeActiveTiles();

        openvdb::tree::IteratorRange<openvdb::BoolGrid::ValueOnCIter> range(border_mask->cbeginValueOn());
        tbb::parallel_for(range, [&](openvdb::tree::IteratorRange<openvdb::BoolGrid::ValueOnCIter> range) {
            auto targetAccessor = resampledLevelSet->getAccessor();
            auto u_weights_accessor = uWeights->getAccessor();
            auto v_weights_accessor = vWeights->getAccessor();
            auto w_weights_accessor = wWeights->getAccessor();

            for (; range; ++range) {
                auto coord = range.iterator().getCoord();
                auto point_000 = targetAccessor.getValue(coord.offsetBy(0, 0, 0));
                auto point_001 = targetAccessor.getValue(coord.offsetBy(0, 0, 1));
                auto point_010 = targetAccessor.getValue(coord.offsetBy(0, 1, 0));
                auto point_011 = targetAccessor.getValue(coord.offsetBy(0, 1, 1));
                auto point_100 = targetAccessor.getValue(coord.offsetBy(1, 0, 0));
                auto point_101 = targetAccessor.getValue(coord.offsetBy(1, 0, 1));
                auto point_110 = targetAccessor.getValue(coord.offsetBy(1, 1, 0));

                u_weights_accessor.setValue(coord, fraction_inside(point_000, point_010, point_001, point_011));
                v_weights_accessor.setValue(coord, fraction_inside(point_000, point_001, point_100, point_101));
                w_weights_accessor.setValue(coord, fraction_inside(point_000, point_010, point_100, point_110));
            }
        });
    };

    computationLambda(domain.fluidLevelSet(), grid.uFluidWeights(), grid.vFluidWeights(), grid.wFluidWeights());
    computationLambda(domain.solidLevelSet(), grid.uWeights(), grid.vWeights(), grid.wWeights());

    // Combine the face are fractions of the solid objects and solid boundary by adding the values
    for (auto solidObj : domain.solidObjs()) {
        openvdb::tools::compSum(*grid.uWeights(), *solidObj->uWeights()->deepCopy());
        openvdb::tools::compSum(*grid.vWeights(), *solidObj->vWeights()->deepCopy());
        openvdb::tools::compSum(*grid.wWeights(), *solidObj->wWeights()->deepCopy());
    }

    // Constrain the face area fraction values
    openvdb::tools::foreach (grid.uWeights()->beginValueAll(),
                             [](const openvdb::FloatGrid::ValueAllIter &iter) { iter.setValue(clamp(iter.getValue(), 0.f, 1.f)); });
    openvdb::tools::foreach (grid.vWeights()->beginValueAll(),
                             [](const openvdb::FloatGrid::ValueAllIter &iter) { iter.setValue(clamp(iter.getValue(), 0.f, 1.f)); });
    openvdb::tools::foreach (grid.wWeights()->beginValueAll(),
                             [](const openvdb::FloatGrid::ValueAllIter &iter) { iter.setValue(clamp(iter.getValue(), 0.f, 1.f)); });

    auto t_end = high_resolution_clock::now();
    auto weights_duration = duration_cast<milliseconds>(t_end - t_start);
}

Eigen::VectorXd FluidSimulator::compute_pressure_divirgence_constraint(
    FluidDomain &domain, openvdb::Int32Grid::Accessor &fluid_indices_accessor, int system_size, float dt) {
    auto &grid = domain.grid();
    // Solver of linear system
    Eigen::setNbThreads(8);
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double, 1 >, Eigen::UpLoType::Lower | Eigen::UpLoType::Upper > solver;
    solver.setTolerance(1e-15);

    // Laplacian matrix
    Eigen::SparseMatrix<double, 1> alpha_matrix;

    alpha_matrix.resize(system_size, system_size);
    alpha_matrix.reserve(Eigen::VectorXi::Constant(system_size, 7));

    for (auto solidObj : domain.solidObjs()) {
        solidObj->jMatrix().resize(6, system_size);
        solidObj->jMatrix().reserve(Eigen::VectorXi::Constant(6, system_size));
    }

    Eigen::VectorXd rhs(system_size);

    auto t_start = high_resolution_clock::now();

    auto solidBbox = domain.solidLevelSet().getActiveCoordBBox();
    auto u_weights_accessor = grid.uWeights()->getAccessor();
    auto v_weights_accessor = grid.vWeights()->getAccessor();
    auto w_weights_accessor = grid.wWeights()->getAccessor();
    auto vel_front_accessor = grid.velFront()->getAccessor();

    auto fluidAccessor = domain.fluidLevelSet().getAccessor();
    auto solidAccessor = domain.solidLevelSet().getAccessor();

    auto mask = openvdb::tools::createPointMask(domain.particleSet(), *domain._i2w_transform);
    openvdb::tools::dilateActiveValues(mask->tree(), 2, openvdb::tools::NN_FACE_EDGE_VERTEX, openvdb::tools::TilePolicy::EXPAND_TILES);

    openvdb::tree::IteratorRange<openvdb::MaskGrid::ValueOnCIter> fluidRange(mask->cbeginValueOn());

    float scale = dt / (FluidDomain::density * domain.voxelSize());

    // Build the system's matrices and right hand side
    for (; fluidRange; ++fluidRange) {
        auto coord = fluidRange.iterator().getCoord();

        if (fluid_indices_accessor.getValue(coord) != -1) {
            int idx = fluid_indices_accessor.getValue(coord);
            rhs[idx] = 0;
            double diagonal = 0;
            double theta = 0;
            double term;

            // These phi values are used to determine in which class the neighbouring cells lie (empty, fluid, solid)
            float phi_i_j_k = fluidAccessor.getValue(coord);
            float phi_i_minus_1_j_k = fluidAccessor.getValue(coord.offsetBy(-1, 0, 0));
            float phi_i_j_minus_1_k = fluidAccessor.getValue(coord.offsetBy(0, -1, 0));
            float phi_i_j_k_minus_1 = fluidAccessor.getValue(coord.offsetBy(0, 0, -1));
            float phi_i_plus_1_j_k = fluidAccessor.getValue(coord.offsetBy(1, 0, 0));
            float phi_i_j_plus_1_k = fluidAccessor.getValue(coord.offsetBy(0, 1, 0));
            float phi_i_j_k_plus_1 = fluidAccessor.getValue(coord.offsetBy(0, 0, 1));

            float solid_phi_i_j_k = solidAccessor.getValue(coord);
            float solid_phi_i_minus_1_j_k = solidAccessor.getValue(coord.offsetBy(-1, 0, 0));
            float solid_phi_i_j_minus_1_k = solidAccessor.getValue(coord.offsetBy(0, -1, 0));
            float solid_phi_i_j_k_minus_1 = solidAccessor.getValue(coord.offsetBy(0, 0, -1));
            float solid_phi_i_plus_1_j_k = solidAccessor.getValue(coord.offsetBy(1, 0, 0));
            float solid_phi_i_j_plus_1_k = solidAccessor.getValue(coord.offsetBy(0, 1, 0));
            float solid_phi_i_j_k_plus_1 = solidAccessor.getValue(coord.offsetBy(0, 0, 1));

            // The fluid face area fraction is 1 - the calculated (solid) face are fraction.

            // Back
            auto weight = (1 - u_weights_accessor.getValue(coord));
            term = weight * scale;
            if (phi_i_minus_1_j_k < 0) {
                alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(-1, 0, 0))) = -term;
                diagonal += term;
            } else if (solid_phi_i_minus_1_j_k > 0) {
                theta = max(0.001f, fraction_inside(phi_i_minus_1_j_k, phi_i_j_k));
                diagonal += term * (1.f / theta);
            }
            rhs[idx] += weight * grid.velHalfIndexed(vel_front_accessor, coord)[0];

            // Front
            weight = (1 - u_weights_accessor.getValue(coord.offsetBy(1, 0, 0)));
            term = weight * scale;
            if (phi_i_plus_1_j_k < 0) {
                alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(1, 0, 0))) = -term;
                diagonal += term;
            } else if (solid_phi_i_plus_1_j_k > 0) {
                theta = max(0.001f, fraction_inside(phi_i_j_k, phi_i_plus_1_j_k));
                diagonal += term * (1.f / theta);
            }
            rhs[idx] -= weight * grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(1, 0, 0))[0];

            // Bot
            weight = (1 - v_weights_accessor.getValue(coord));
            term = weight * scale;
            if (phi_i_j_minus_1_k < 0) {
                alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(0, -1, 0))) = -term;
                diagonal += term;
            } else if (solid_phi_i_j_minus_1_k > 0) {
                theta = max(0.001f, fraction_inside(phi_i_j_minus_1_k, phi_i_j_k));
                diagonal += term * (1.f / theta);
            }
            rhs[idx] += weight * grid.velHalfIndexed(vel_front_accessor, coord)[1];

            // Top
            weight = (1 - v_weights_accessor.getValue(coord.offsetBy(0, 1, 0)));
            term = weight * scale;
            if (phi_i_j_plus_1_k < 0) {
                alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(0, 1, 0))) = -term;
                diagonal += term;
            } else if (solid_phi_i_j_plus_1_k > 0) {
                theta = max(0.001f, fraction_inside(phi_i_j_k, phi_i_j_plus_1_k));
                diagonal += term * (1.f / theta);
            }
            rhs[idx] -= weight * grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 1, 0))[1];

            // Left
            weight = (1 - w_weights_accessor.getValue(coord));
            term = weight * scale;
            if (phi_i_j_k_minus_1 < 0) {
                alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(0, 0, -1))) = -term;
                diagonal += term;
            } else if (solid_phi_i_j_k_minus_1 > 0) {
                theta = max(0.001f, fraction_inside(phi_i_j_k_minus_1, phi_i_j_k));
                diagonal += term * (1.f / theta);
            }
            rhs[idx] += weight * grid.velHalfIndexed(vel_front_accessor, coord)[2];

            // Right
            weight = (1 - w_weights_accessor.getValue(coord.offsetBy(0, 0, 1)));
            term = weight * scale;
            if (phi_i_j_k_plus_1 < 0) {
                alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(0, 0, 1))) = -term;
                diagonal += term;
            } else if (solid_phi_i_j_k_plus_1 > 0) {
                theta = max(0.001f, fraction_inside(phi_i_j_k, phi_i_j_k_plus_1));
                diagonal += term * (1.f / theta);
            }
            rhs[idx] -= weight * grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 0, 1))[2];

            alpha_matrix.insert(idx, idx) = diagonal;

            // J matrix building for strong two-way (solid-fluid) coupling
            for (auto solidObj : domain.solidObjs()) {
                solidObj->constructJMatrix(idx, coord, domain.voxelSize());
            }
        }
    }
    // });
    alpha_matrix *= sqr(domain.voxelSize());
    rhs *= sqr(domain.voxelSize());

    // Integrate the contribution of solid objects in the alpha matrix and rhs
    // for strong solid-fluid coupling pressure solve
    for (auto solidObj : domain.solidObjs()) {
        auto &jMatrix = solidObj->jMatrix();
        jMatrix.prune(0.);
        alpha_matrix += dt * jMatrix.transpose() * solidObj->massMatrixInverse() * jMatrix;

        rhs -= jMatrix.transpose() * solidObj->combinedVelocityVector();
    }

    auto t_end = high_resolution_clock::now();
    auto system_build_duration = duration_cast<milliseconds>(t_end - t_start);

    Eigen::VectorXd pressure_grid(system_size);
    solver.compute(alpha_matrix);
    pressure_grid = solver.solve(rhs);
    #ifdef PRINT
    std::cout << rhs << std::endl;
    std::cout << pressure_grid << std::endl;
    #endif
    std::cout << "\n" << solver.error() << " " << solver.info() << " " << solver.iterations() << "\n" << std::endl;
    return pressure_grid;
}
    auto system_build_duration = duration_cast<milliseconds>(t_end - t_start);

    t_start = high_resolution_clock::now();
    Eigen::VectorXd pressure_grid(system_size);
    solver.compute(alpha_matrix);
    pressure_grid = solver.solve(rhs);
    // std::cout << "\n" << solver.error() << " " << solver.info() << " " << solver.iterations() << "\n" << std::endl;
    t_end = high_resolution_clock::now();
    auto solver_duration = duration_cast<milliseconds>(t_end - t_start);
    // std::cout << std::endl;
    // for (auto it = rhs.begin(); it != rhs.end(); ++it) { std::cout << *it << " "; }
    // std::cout << std::endl;
    // std::cout << std::endl;
    // for (auto it = pressure_grid.begin(); it != pressure_grid.end(); ++it) { std::cout << *it << " "; }
    // std::cout << std::endl;

void FluidSimulator::project_divirgence_constraint(FluidDomain &domain,
                                                   openvdb::Int32Grid::Accessor &fluid_indices_accessor,
                                                   Eigen::VectorXd &pressure_grid, float dt) {
    auto &grid = domain.grid();

    // Project velocities based on the calculated pressure
    domain.grid().validUFront()->clear();
    domain.grid().validVFront()->clear();
    domain.grid().validWFront()->clear();

    auto u_valid_acc = domain.grid().validUFront()->getAccessor();
    auto v_valid_acc = domain.grid().validVFront()->getAccessor();
    auto w_valid_acc = domain.grid().validWFront()->getAccessor();

    auto fluidAccessor = domain.fluidLevelSet().getAccessor();
    auto solidAccessor = domain.solidLevelSet().getAccessor();

    auto u_weights_accessor = grid.uWeights()->getAccessor();
    auto v_weights_accessor = grid.vWeights()->getAccessor();
    auto w_weights_accessor = grid.wWeights()->getAccessor();
    auto vel_front_accessor = grid.velFront()->getAccessor();
    auto vel_back_accessor = grid.velBack()->getAccessor();

    auto solidBbox = domain.solidLevelSet().getActiveCoordBBox();

    grid.velBack()->clear();

    vel_back_accessor = grid.velBack()->getAccessor();
    for (auto it = solidBbox.beginZYX(); it != solidBbox.endZYX(); ++it) {
        auto coord = (*it);

        int idx = fluid_indices_accessor.getValue(coord);
        int idx_i_minus1 = fluid_indices_accessor.getValue(coord.offsetBy(-1, 0, 0));
        int idx_j_minus1 = fluid_indices_accessor.getValue(coord.offsetBy(0, -1, 0));
        int idx_k_minus1 = fluid_indices_accessor.getValue(coord.offsetBy(0, 0, -1));

        float p = idx >= 0 ? pressure_grid[idx] : 0;
        float p_i_minus1 = idx_i_minus1 >= 0 ? pressure_grid[idx_i_minus1] : 0;
        float p_j_minus1 = idx_j_minus1 >= 0 ? pressure_grid[idx_j_minus1] : 0;
        float p_k_minus1 = idx_k_minus1 >= 0 ? pressure_grid[idx_k_minus1] : 0;

        float phi_i_j_k = fluidAccessor.getValue(coord);
        float phi_i_minus_1_j_k = fluidAccessor.getValue(coord.offsetBy(-1, 0, 0));
        float phi_i_j_minus_1_k = fluidAccessor.getValue(coord.offsetBy(0, -1, 0));
        float phi_i_j_k_minus_1 = fluidAccessor.getValue(coord.offsetBy(0, 0, -1));

        auto prevVel = grid.velHalfIndexed(vel_front_accessor, coord);
        openvdb::Vec3d newVel = prevVel;

        if ((1 - u_weights_accessor.getValue(coord)) > 0 && (phi_i_j_k < 0 || phi_i_minus_1_j_k < 0)) {
            float theta = 1;
            if (phi_i_j_k >= 0 || phi_i_minus_1_j_k >= 0) {
                theta = max(0.001f, fraction_inside(phi_i_minus_1_j_k, phi_i_j_k));
            }
            newVel[0] -= dt / FluidDomain::density * (p - p_i_minus1) / domain.voxelSize() / theta;
            u_valid_acc.setValueOn(coord);
        }

        if ((1 - v_weights_accessor.getValue(coord)) > 0 && (phi_i_j_k < 0 || phi_i_j_minus_1_k < 0)) {
            float theta = 1;
            if (phi_i_j_k >= 0 || phi_i_j_minus_1_k >= 0) {
                theta = max(0.001f, fraction_inside(phi_i_j_minus_1_k, phi_i_j_k));
            }
            newVel[1] -= dt / FluidDomain::density * (p - p_j_minus1) / domain.voxelSize() / theta;
            v_valid_acc.setValueOn(coord);
        }

        if ((1 - w_weights_accessor.getValue(coord)) > 0 && (phi_i_j_k < 0 || phi_i_j_k_minus_1 < 0)) {
            float theta = 1;
            if (phi_i_j_k >= 0 || phi_i_j_k_minus_1 >= 0) {
                theta = max(0.001f, fraction_inside(phi_i_j_k_minus_1, phi_i_j_k));
            }
            newVel[2] -= dt / FluidDomain::density * (p - p_k_minus1) / domain.voxelSize() / theta;
            w_valid_acc.setValueOn(coord);
        }
        grid.setVelHalfIndexed(vel_back_accessor, coord, newVel);
    }

    for (auto solidObj : domain.solidObjs()) {
        auto new_solid_vel_vector = solidObj->combinedVelocityVector()
                                    + dt * solidObj->massMatrixInverse() * solidObj->jMatrix() * pressure_grid;
        // Eigen::Vector<double, 6> force = solidObj->massMatrixInverse() * solidObj->jMatrix() * pressure_grid;
        solidObj->setVelocityVector(new_solid_vel_vector);
    }
    grid.swapVelocityBuffers();
}

    t_end = high_resolution_clock::now();
    auto update_duration = duration_cast<milliseconds>(t_end - t_start);

    printf("weights_duration %.3f, system_build_duration %.3f, solver_duration %.3f, update_duration %.3f \n",
           weights_duration.count() / 1000.f, system_build_duration.count() / 1000.f, solver_duration.count() /
           1000.f, update_duration.count() / 1000.f);

    grid.swapVelocityBuffers();
}

void FluidSimulator::constrain_velocity(FluidDomain &domain) {
    auto &grid = domain.grid();
    auto gradient = openvdb::tools::gradient(*domain.solidLevelSet().getLevelSet());
    grid.setVelBack(grid.velFront()->deepCopy());

    grid.uWeights()->tree().topologyIntersection(domain.solidLevelSet().getLevelSet()->tree());
    grid.vWeights()->tree().topologyIntersection(domain.solidLevelSet().getLevelSet()->tree());
    grid.wWeights()->tree().topologyIntersection(domain.solidLevelSet().getLevelSet()->tree());

    auto u_weights_accessor = grid.uWeights()->getAccessor();
    auto v_weights_accessor = grid.vWeights()->getAccessor();
    auto w_weights_accessor = grid.wWeights()->getAccessor();
    auto vel_back_accessor = grid.velBack()->getAccessor();

    openvdb::tools::GridSampler<openvdb::Vec3SGrid::Accessor, openvdb::tools::BoxSampler> grad_staggered_x_sampler(
        gradient->getAccessor(), gradient->transform());
    openvdb::tools::GridSampler<openvdb::Vec3SGrid::Accessor, openvdb::tools::BoxSampler> grad_staggered_y_sampler(
        gradient->getAccessor(), gradient->transform());
    openvdb::tools::GridSampler<openvdb::Vec3SGrid::Accessor, openvdb::tools::BoxSampler> grad_staggered_z_sampler(
        gradient->getAccessor(), gradient->transform());

    MacGrid::StaggeredSampler vel_staggered_x_sampler(grid.velFront()->getAccessor(), grid.velFront()->transform());
    MacGrid::StaggeredSampler vel_staggered_y_sampler(grid.velFront()->getAccessor(), grid.velFront()->transform());
    MacGrid::StaggeredSampler vel_staggered_z_sampler(grid.velFront()->getAccessor(), grid.velFront()->transform());

    auto unionMask = openvdb::createGrid<openvdb::MaskGrid>();
    unionMask->tree().topologyUnion(grid.uWeights()->tree());
    unionMask->tree().topologyUnion(grid.vWeights()->tree());
    unionMask->tree().topologyUnion(grid.wWeights()->tree());

    for (auto iter = unionMask->beginValueOn(); iter; ++iter) {
        auto coord = iter.getCoord();
        openvdb::Vec3d newVel = vel_back_accessor.getValue(coord);
        if (u_weights_accessor.getValue(coord) == 1) {
            // auto normal = grad_staggered_x_accessor.getValue(coord).unit();
            auto normal = grad_staggered_x_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0.5, 0, 0));
            if (!normal.isZero()) normal.normalize();
            // auto vel = vel_staggered_x_accessor.getValue(coord);
            auto vel = vel_staggered_x_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0.5, 0, 0));
            auto projection = vel.projection(normal);
            auto tangential = vel - projection;
            newVel[0] = tangential[0];
        }
        if (v_weights_accessor.getValue(coord) == 1) {
            // auto normal = grad_staggered_y_accessor.getValue(coord).unit();
            auto normal = grad_staggered_y_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0.5, 0));
            if (!normal.isZero()) normal.normalize();
            // auto vel = vel_staggered_y_accessor.getValue(coord);
            auto vel = vel_staggered_y_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0.5, 0));
            auto projection = vel.projection(normal);
            auto tangential = vel - projection;
            newVel[1] = tangential[1];
        }
        if (w_weights_accessor.getValue(coord) == 1) {
            // auto normal = grad_staggered_z_accessor.getValue(coord).unit();
            auto normal = grad_staggered_z_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0, 0.5));
            if (!normal.isZero()) normal.normalize();
            // auto vel = vel_staggered_z_accessor.getValue(coord);
            auto vel = vel_staggered_z_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0, 0.5));
            auto projection = vel.projection(normal);
            auto tangential = vel - projection;
            newVel[2] = tangential[2];
        }
        grid.setVelHalfIndexed(vel_back_accessor, coord, newVel);
    }

    grid.swapVelocityBuffers();
}

void FluidSimulator::constrain_velocity(FluidDomain &domain, SolidObject solidObject) {
    auto &grid = domain.grid();

    auto vel_trans = solidObject.velTranslational();
    auto vel_ang = solidObject.velAngular();

    auto gradient = openvdb::tools::gradient(*solidObject.levelSet().getLevelSet());
    grid.setVelBack(grid.velFront()->deepCopy());

    solidObject.uWeights()->tree().topologyIntersection(solidObject.levelSet().getLevelSet()->tree());
    solidObject.vWeights()->tree().topologyIntersection(solidObject.levelSet().getLevelSet()->tree());
    solidObject.wWeights()->tree().topologyIntersection(solidObject.levelSet().getLevelSet()->tree());

    auto u_weights_accessor = solidObject.uWeights()->getAccessor();
    auto v_weights_accessor = solidObject.vWeights()->getAccessor();
    auto w_weights_accessor = solidObject.wWeights()->getAccessor();
    auto vel_back_accessor = grid.velBack()->getAccessor();

    openvdb::tools::GridSampler<openvdb::Vec3SGrid::Accessor, openvdb::tools::BoxSampler> grad_staggered_x_sampler(
        gradient->getAccessor(), gradient->transform());
    openvdb::tools::GridSampler<openvdb::Vec3SGrid::Accessor, openvdb::tools::BoxSampler> grad_staggered_y_sampler(
        gradient->getAccessor(), gradient->transform());
    openvdb::tools::GridSampler<openvdb::Vec3SGrid::Accessor, openvdb::tools::BoxSampler> grad_staggered_z_sampler(
        gradient->getAccessor(), gradient->transform());

    MacGrid::StaggeredSampler vel_staggered_x_sampler(grid.velFront()->getAccessor(), grid.velFront()->transform());
    MacGrid::StaggeredSampler vel_staggered_y_sampler(grid.velFront()->getAccessor(), grid.velFront()->transform());
    MacGrid::StaggeredSampler vel_staggered_z_sampler(grid.velFront()->getAccessor(), grid.velFront()->transform());

    auto unionMask = openvdb::createGrid<openvdb::MaskGrid>();
    unionMask->tree().topologyUnion(solidObject.uWeights()->tree());
    unionMask->tree().topologyUnion(solidObject.vWeights()->tree());
    unionMask->tree().topologyUnion(solidObject.wWeights()->tree());

    for (auto iter = unionMask->beginValueOn(); iter; ++iter) {
        auto coord = iter.getCoord();
        openvdb::Vec3d newVel = vel_back_accessor.getValue(coord);
        if (u_weights_accessor.getValue(coord) == 1) {
            auto solidVel = vel_trans + domain.voxelSize() * vel_ang.cross(coord.asVec3d() - openvdb::Vec3d(0.5, 0, 0));
            // auto normal = grad_staggered_x_accessor.getValue(coord).unit();
            auto normal = grad_staggered_x_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0.5, 0, 0));
            if (!normal.isZero()) normal.normalize();
            auto vel = vel_staggered_x_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0.5, 0, 0));
            auto projection = vel.projection(normal);
            auto tangential = vel - projection;
            newVel[0] = tangential[0] + normal.dot(solidVel) * normal[0];
        }
        if (v_weights_accessor.getValue(coord) == 1) {
            auto solidVel = vel_trans + domain.voxelSize() * vel_ang.cross(coord.asVec3d() - openvdb::Vec3d(0, 0.5, 0));
            // auto normal = grad_staggered_y_accessor.getValue(coord).unit();
            auto normal = grad_staggered_y_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0.5, 0));
            if (!normal.isZero()) normal.normalize();
            auto vel = vel_staggered_y_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0.5, 0));
            auto projection = vel.projection(normal);
            auto tangential = vel - projection;
            newVel[1] = tangential[1] + normal.dot(solidVel) * normal[1];
        }
        if (w_weights_accessor.getValue(coord) == 1) {
            auto solidVel = vel_trans + domain.voxelSize() * vel_ang.cross(coord.asVec3d() - openvdb::Vec3d(0, 0, 0.5));

            // auto normal = grad_staggered_z_accessor.getValue(coord).unit();
            auto normal = grad_staggered_z_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0, 0.5));
            if (!normal.isZero()) normal.normalize();
            auto vel = vel_staggered_z_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0, 0.5));
            auto projection = vel.projection(normal);
            auto tangential = vel - projection;
            newVel[2] = tangential[2] + normal.dot(solidVel) * normal[2];
        }
        grid.setVelHalfIndexed(vel_back_accessor, coord, newVel);
    }

    grid.swapVelocityBuffers();
}

double FluidSimulator::calculate_kernel_function(openvdb::Vec3d vec) {
    double val1 = calculate_trilinear_hat(vec[0]);
    double val2 = calculate_trilinear_hat(vec[1]);
    double val3 = calculate_trilinear_hat(vec[2]);
    return val1 * val2 * val3;
}

double FluidSimulator::calculate_kernel_function(double x, double y, double z) {
    double val1 = calculate_trilinear_hat(x);
    double val2 = calculate_trilinear_hat(y);
    double val3 = calculate_trilinear_hat(z);
    return val1 * val2 * val3;
}

openvdb::Vec3d FluidSimulator::calculate_kernel_function_staggered(openvdb::Vec3d difference) {
    float val_x_0 = calculate_trilinear_hat(difference[0]);
    float val_y_0 = calculate_trilinear_hat(difference[1]);
    float val_z_0 = calculate_trilinear_hat(difference[2]);

    float val_x_1 = calculate_trilinear_hat(difference[0] + 0.5);
    float val_y_1 = calculate_trilinear_hat(difference[1] + 0.5);
    float val_z_1 = calculate_trilinear_hat(difference[2] + 0.5);

    return openvdb::Vec3d(val_x_1 * val_y_0 * val_z_0, val_x_0 * val_y_1 * val_z_0, val_x_0 * val_y_0 * val_z_1);
}