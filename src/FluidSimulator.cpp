#pragma once
#include "FluidSimulator.h"

#include "FluidDomain.h"
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

class GatherTransfer {
public:
    const openvdb::Vec3dGrid::Ptr &vel_weight;
    const openvdb::Vec3dGrid::Ptr &vel;
    openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &p_atlas;
    openvdb::CoordBBox &bbox;
    FluidDomain &domain;
    void operator()(const tbb::blocked_range<int> &i_range) const {
        auto back_accessor = vel->getAccessor();
        auto vel_weight_accessor = vel_weight->getAccessor();
        openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Iterator p_it(*p_atlas);
        for (int i = i_range.begin(); i < i_range.end(); ++i) {
            for (int j = bbox.min()[1]; j < bbox.max()[1]; ++j) {
                for (int k = bbox.min()[2]; k < bbox.max()[2]; ++k) {
                    p_it.worldSpaceSearchAndUpdate(
                        domain.i2w_transform->indexToWorld(
                            openvdb::BBoxd(openvdb::Vec3d(i - 1, j - 1, k - 1), openvdb::Vec3d(i + 1, j + 1, k + 1))),
                        domain.particleSet());
                    for (; p_it; ++p_it) {
                        auto &p = domain.particleSet()[*p_it];

                        auto coord = openvdb::Coord(i, j, k);
                        openvdb::Vec3d weight_vec(
                            FluidSimulator::calculate_kernel_function_staggered(p.pos() - coord.asVec3d()));
                        if (!weight_vec.isZero()) {
                            vel_weight_accessor.setValue(coord, vel_weight_accessor.getValue(coord) + weight_vec);
                            back_accessor.setValue(coord, back_accessor.getValue(coord) + weight_vec * p.vel());
                        }
                    }
                }
            }
        }
    }

    GatherTransfer(FluidDomain &domain, openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &p_atlas,
                   openvdb::CoordBBox &bbox, const openvdb::Vec3dGrid::Ptr &vel,
                   const openvdb::Vec3dGrid::Ptr &vel_weight):
     domain(domain),
     p_atlas(p_atlas), bbox(bbox), vel_weight(vel_weight), vel(vel) {}
};

class ShootingTransfer {
public:
    openvdb::Vec3dGrid::Ptr vel_weight;
    openvdb::Vec3dGrid::Ptr vel;
    ParticleSet &p_set;
    FluidDomain &domain;
    void operator()(const tbb::blocked_range<int> &p_range) {
        auto back_accessor = vel->getAccessor();
        auto vel_weight_accessor = vel_weight->getAccessor();
        for (int p = p_range.begin(); p < p_range.end(); ++p) {
            auto &particle = p_set[p];
            auto start_coord = openvdb::Coord::round(particle.pos()).offsetBy(-1, -1, -1);
            auto end_coord = start_coord.offsetBy(3, 3, 3);
            openvdb::CoordBBox bbox(start_coord, end_coord);

            for (int i = bbox.min()[0]; i < bbox.max()[0]; ++i) {
                for (int j = bbox.min()[1]; j < bbox.max()[1]; ++j) {
                    for (int k = bbox.min()[2]; k < bbox.max()[2]; ++k) {
                        auto coord = openvdb::Coord(i, j, k);
                        openvdb::Vec3d weight_vec(
                            FluidSimulator::calculate_kernel_function_staggered(particle.pos() - coord.asVec3d()));
                        if (!weight_vec.isZero()) {
                            vel_weight_accessor.setValue(coord, vel_weight_accessor.getValue(coord) + weight_vec);
                            back_accessor.setValue(coord, back_accessor.getValue(coord) + weight_vec * particle.vel());
                        }
                    }
                }
            }
        }
    }

    ShootingTransfer(ShootingTransfer &x, tbb::split):
     domain(x.domain), p_set(x.p_set), vel_weight(openvdb::Vec3DGrid::create(openvdb::Vec3d(epsilon))),
     vel(openvdb::Vec3DGrid::create(openvdb::Vec3d(0))) {}

    ShootingTransfer(FluidDomain &domain, ParticleSet &p_set):
     domain(domain), p_set(p_set), vel_weight(openvdb::Vec3DGrid::create(openvdb::Vec3d(epsilon))),
     vel(openvdb::Vec3DGrid::create(openvdb::Vec3d(0))) {
        vel_weight->setTransform(domain.i2w_transform);
        vel->setGridClass(openvdb::GRID_STAGGERED);
        vel->setTransform(domain.i2w_transform);
    }
    void join(const ShootingTransfer &y) {
        openvdb::tools::compSum(*vel, *y.vel);
        openvdb::tools::compSum(*vel_weight, *y.vel_weight);
    }
};

class ReseedingFunctor {
    std::random_device rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen;      // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis;
    const int max_particles;
    const int min_particles;
    openvdb::MaskGrid::Ptr active_mask;
    openvdb::Vec3dGrid::Ptr vel;
    FluidDomain &domain;

public:
    std::vector<int> particles_to_be_deleted;
    std::vector<Particle> particles_to_be_added;

    openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &p_atlas;
    void operator()(openvdb::tree::IteratorRange<openvdb::MaskGrid::ValueOnCIter> &range) {
        MacGrid::BoxSampler vel_sampler(vel->getAccessor(), vel->transform());
        openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Iterator p_it(*p_atlas);
        auto active_mask_accessor = active_mask->getAccessor();
        for (; range; ++range) {
            auto &it = range.iterator();
            auto cell_bbox = openvdb::BBoxd(it.getCoord().asVec3d() - openvdb::Vec3d(0.5, 0.5, 0.5),
                                            it.getCoord().asVec3d() + openvdb::Vec3d(0.5, 0.5, 0.5));

            p_it.worldSpaceSearchAndUpdate(domain.i2w_transform->indexToWorld(cell_bbox), domain.particleSet());

            openvdb::BBoxd particles_bbox = openvdb::BBoxd(openvdb::Vec3d(1, 1, 1), openvdb::Vec3d(-1, -1, -1));
            std::vector<int> particles_inside;
            int p_counter = 0;
            int test_counter = 0;
            for (; p_it; ++p_it) {
                auto &p = domain.particleSet()[*p_it];
                test_counter++;

                if (cell_bbox.isInside(p.pos())) {
                    particles_inside.push_back(*p_it);
                    p_counter++;
                    // Initialise bbox
                    if (particles_bbox.max() < particles_bbox.min()) {
                        particles_bbox = openvdb::BBoxd(p.pos(), p.pos());
                    } else {
                        particles_bbox.expand(p.pos());
                    }
                }
            }
            if (p_counter > max_particles) {
                for (auto it = particles_inside.begin(); p_counter > max_particles && it != particles_inside.end();
                     ++it) {
                    auto particles_bbox_min = particles_bbox.min();
                    auto particles_bbox_max = particles_bbox.max();
                    auto particle_pos = domain.particleSet()[*it].pos();
                    if (particles_bbox_min[0] < particle_pos[0] && particles_bbox_min[1] < particle_pos[1]
                        && particles_bbox_min[2] < particle_pos[2] && particle_pos[0] < particles_bbox_max[0]
                        && particle_pos[1] < particles_bbox_max[1] && particle_pos[2] < particles_bbox_max[2]) {
                        particles_to_be_deleted.push_back(*it);
                        p_counter--;
                    }
                }
            } else if (p_counter < min_particles) {
                if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(-1, 0, 0))) {
                    auto new_extent = cell_bbox.min();
                    new_extent[0] = particles_bbox.min()[0];
                    cell_bbox = openvdb::BBoxd(new_extent, cell_bbox.max());
                }
                if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(0, -1, 0))) {
                    auto new_extent = cell_bbox.min();
                    new_extent[1] = particles_bbox.min()[1];
                    cell_bbox = openvdb::BBoxd(new_extent, cell_bbox.max());
                }
                if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(0, 0, -1))) {
                    auto new_extent = cell_bbox.min();
                    new_extent[2] = particles_bbox.min()[2];
                    cell_bbox = openvdb::BBoxd(new_extent, cell_bbox.max());
                }
                if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(1, 0, 0))) {
                    auto new_extent = cell_bbox.max();
                    new_extent[0] = particles_bbox.max()[0];
                    cell_bbox = openvdb::BBoxd(cell_bbox.min(), new_extent);
                }
                if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(0, 1, 0))) {
                    auto new_extent = cell_bbox.max();
                    new_extent[1] = particles_bbox.max()[1];
                    cell_bbox = openvdb::BBoxd(cell_bbox.min(), new_extent);
                }
                if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(0, 0, 1))) {
                    auto new_extent = cell_bbox.max();
                    new_extent[2] = particles_bbox.max()[2];
                    cell_bbox = openvdb::BBoxd(cell_bbox.min(), new_extent);
                }

                while (p_counter++ < min_particles) {
                    openvdb::Vec3d new_pos =
                        cell_bbox.min()
                        + openvdb::Vec3d(dis(gen), dis(gen), dis(gen)) * (cell_bbox.max() - cell_bbox.min());
                    openvdb::Vec3d new_vel = domain.grid().velInterpolatedI(vel_sampler, new_pos);
                    Particle p(new_pos, new_vel);
                    particles_to_be_added.push_back(p);
                }
            }
        }
    }

    ReseedingFunctor(ReseedingFunctor &x, tbb::split):
     domain(x.domain), p_atlas(x.p_atlas), active_mask(x.active_mask), vel(x.vel), min_particles(x.min_particles),
     max_particles(x.max_particles) {
        gen = std::mt19937(rd());
        dis = std::uniform_real_distribution(0., 1.);
    }

    ReseedingFunctor(FluidDomain &domain, openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &p_atlas,
                     openvdb::MaskGrid::Ptr active_mask, openvdb::Vec3dGrid::Ptr vel, const int min_particles,
                     const int max_particles):
     domain(domain),
     p_atlas(p_atlas), active_mask(active_mask), vel(vel), min_particles(min_particles), max_particles(max_particles) {
        gen = std::mt19937(rd());
        dis = std::uniform_real_distribution(0., 1.);
    }

    void join(const ReseedingFunctor &y) {
        particles_to_be_added.reserve(particles_to_be_added.size() + y.particles_to_be_added.size());
        particles_to_be_deleted.reserve(particles_to_be_deleted.size() + y.particles_to_be_deleted.size());

        particles_to_be_added.insert(particles_to_be_added.end(), y.particles_to_be_added.begin(),
                                     y.particles_to_be_added.end());
        particles_to_be_deleted.insert(particles_to_be_deleted.end(), y.particles_to_be_deleted.begin(),
                                       y.particles_to_be_deleted.end());
    }
};

FluidSimulator::FluidSimulator(openvdb::CoordBBox bbox): _bbox(bbox) {}

FluidSimulator::~FluidSimulator() {}

void FluidSimulator::print_velocity_field(MacGrid &grid, const char *variable_name) {
#ifdef PRINT
    auto accessor = grid.velFront()->getAccessor();
    auto bbox = grid.velFront()->evalActiveVoxelBoundingBox();
    std::cout << variable_name;
    // assert(bbox.min() < bbox.max());
    for (int i = _bbox.min()[0]; i < _bbox.max()[0]; ++i) {
        std::cout << "\n" << i << std::endl;
        for (int j = _bbox.max()[1]; j >= _bbox.min()[1]; --j) {
            std::cout << std::endl;
            for (int k = _bbox.min()[2]; k < _bbox.max()[2]; ++k) {
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
    vel_weight->setTransform(domain.i2w_transform);

    // Zero out radius for the atlas creation as we care only for the particles whose center lies inside a voxel
    auto radius = domain.particleSet().getRadius();
    domain.particleSet().setRadius(0);

    auto p_atlas = openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::create<ParticleSet>(
        domain.particleSet(), domain.voxelSize());

    auto pointMask = openvdb::tools::createPointMask(domain.particleSet(), grid.velBack()->transform());
    openvdb::tools::dilateActiveValues(pointMask->tree(), 2, openvdb::tools::NN_FACE_EDGE_VERTEX);

    grid.velBack()->tree().topologyUnion(pointMask->tree());
    grid.velBack()->tree().voxelizeActiveTiles();

    vel_weight->tree().topologyUnion(pointMask->tree());
    vel_weight->tree().voxelizeActiveTiles();

    auto bbox = pointMask->evalActiveVoxelBoundingBox();
    tbb::parallel_for(tbb::blocked_range<int>(bbox.min()[0], bbox.max()[0], 1),
                      GatherTransfer(domain, p_atlas, bbox, grid.velBack(), vel_weight));
    domain.particleSet().setRadius(radius);
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
        openvdb::tools::GridSampler<openvdb::Vec3dGrid::Accessor, openvdb::tools::StaggeredBoxSampler> vel_sampler(
            domain.grid().velFront()->getAccessor(), domain.grid().velFront()->transform());
        openvdb::tools::GridSampler<openvdb::Vec3dGrid::Accessor, openvdb::tools::StaggeredBoxSampler> vel_diff_sampler(
            domain.grid().velDiff()->getAccessor(), domain.grid().velDiff()->transform());
        for (int i = range.begin(); i < range.end(); ++i) {
            auto &p = domain.particleSet()[i];
            openvdb::Vec3d pic_vel = vel_sampler.isSample(p.pos());
            openvdb::Vec3d flip_vel = p.vel() + vel_diff_sampler.isSample(p.pos());

            p.setVelocity(pic_vel * (1 - flip_pic_ratio) + flip_vel * flip_pic_ratio);
        }
    });
}

void FluidSimulator::reseeding(FluidDomain &domain) {
    float temp_radius = domain.particleSet().getRadius();
    domain.particleSet().setRadius(0);
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
    domain.particleSet().setRadius(temp_radius);
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

        for (auto it = _bbox.beginZYX(); it != _bbox.endZYX(); ++it) {
            auto coord = (*it);

            // Check if this cell will be updated in the x dimension
            if (!valid_mask_u_front_buffer_accessor.isValueOn(coord)) {
                float new_vel_x = 0;
                int n_valid_neighbors_x = 0;

                // Get values of all valid neighbors
                if (valid_mask_u_front_buffer_accessor.isValueOn(coord.offsetBy(-1, 0, 0))) {
                    new_vel_x += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(-1, 0, 0))[0];
                    n_valid_neighbors_x++;
                }
                if (valid_mask_u_front_buffer_accessor.isValueOn(coord.offsetBy(0, -1, 0))) {
                    new_vel_x += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, -1, 0))[0];
                    n_valid_neighbors_x++;
                }
                if (valid_mask_u_front_buffer_accessor.isValueOn(coord.offsetBy(0, 0, -1))) {
                    new_vel_x += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 0, -1))[0];
                    n_valid_neighbors_x++;
                }
                if (valid_mask_u_front_buffer_accessor.isValueOn(coord.offsetBy(1, 0, 0))) {
                    new_vel_x += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(1, 0, 0))[0];
                    n_valid_neighbors_x++;
                }
                if (valid_mask_u_front_buffer_accessor.isValueOn(coord.offsetBy(0, 1, 0))) {
                    new_vel_x += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 1, 0))[0];
                    n_valid_neighbors_x++;
                }
                if (valid_mask_u_front_buffer_accessor.isValueOn(coord.offsetBy(0, 0, 1))) {
                    new_vel_x += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 0, 1))[0];
                    n_valid_neighbors_x++;
                }

                // Average the value for the current cell
                if (n_valid_neighbors_x > 0) {
                    grid.setVelXHalfIndexed(vel_back_accessor, coord, new_vel_x / n_valid_neighbors_x);
                    valid_mask_u_back_buffer_accessor.setValueOn(coord);
                }
            }

            // Check if this cell will be updated in the y dimension
            if (!valid_mask_v_front_buffer_accessor.isValueOn(coord)) {
                float new_vel_y = 0;
                int n_valid_neighbors_y = 0;

                // Get values of all valid neighbors
                if (valid_mask_v_front_buffer_accessor.isValueOn(coord.offsetBy(-1, 0, 0))) {
                    new_vel_y += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(-1, 0, 0))[1];
                    n_valid_neighbors_y++;
                }
                if (valid_mask_v_front_buffer_accessor.isValueOn(coord.offsetBy(0, -1, 0))) {
                    new_vel_y += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, -1, 0))[1];
                    n_valid_neighbors_y++;
                }
                if (valid_mask_v_front_buffer_accessor.isValueOn(coord.offsetBy(0, 0, -1))) {
                    new_vel_y += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 0, -1))[1];
                    n_valid_neighbors_y++;
                }
                if (valid_mask_v_front_buffer_accessor.isValueOn(coord.offsetBy(1, 0, 0))) {
                    new_vel_y += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(1, 0, 0))[1];
                    n_valid_neighbors_y++;
                }
                if (valid_mask_v_front_buffer_accessor.isValueOn(coord.offsetBy(0, 1, 0))) {
                    new_vel_y += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 1, 0))[1];
                    n_valid_neighbors_y++;
                }
                if (valid_mask_v_front_buffer_accessor.isValueOn(coord.offsetBy(0, 0, 1))) {
                    new_vel_y += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 0, 1))[1];
                    n_valid_neighbors_y++;
                }

                // Average the value for the current cell
                if (n_valid_neighbors_y > 0) {
                    grid.setVelYHalfIndexed(vel_back_accessor, coord, new_vel_y / n_valid_neighbors_y);
                    valid_mask_v_back_buffer_accessor.setValueOn(coord);
                }
            }

            if (!valid_mask_w_front_buffer_accessor.isValueOn(coord)) {
                float new_vel_z = 0;
                int n_valid_neighbors_z = 0;

                // Get values of all valid neighbors
                if (valid_mask_w_front_buffer_accessor.isValueOn(coord.offsetBy(-1, 0, 0))) {
                    new_vel_z += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(-1, 0, 0))[2];
                    n_valid_neighbors_z++;
                }
                if (valid_mask_w_front_buffer_accessor.isValueOn(coord.offsetBy(0, -1, 0))) {
                    new_vel_z += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, -1, 0))[2];
                    n_valid_neighbors_z++;
                }
                if (valid_mask_w_front_buffer_accessor.isValueOn(coord.offsetBy(0, 0, -1))) {
                    new_vel_z += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 0, -1))[2];
                    n_valid_neighbors_z++;
                }
                if (valid_mask_w_front_buffer_accessor.isValueOn(coord.offsetBy(1, 0, 0))) {
                    new_vel_z += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(1, 0, 0))[2];
                    n_valid_neighbors_z++;
                }
                if (valid_mask_w_front_buffer_accessor.isValueOn(coord.offsetBy(0, 1, 0))) {
                    new_vel_z += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 1, 0))[2];
                    n_valid_neighbors_z++;
                }
                if (valid_mask_w_front_buffer_accessor.isValueOn(coord.offsetBy(0, 0, 1))) {
                    new_vel_z += grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 0, 1))[2];
                    n_valid_neighbors_z++;
                }

                // Average the value for the current cell
                if (n_valid_neighbors_z > 0) {
                    grid.setVelZHalfIndexed(vel_back_accessor, coord, new_vel_z / n_valid_neighbors_z);
                    valid_mask_w_back_buffer_accessor.setValueOn(coord);
                }
            }
        }
        // Swap valid mask
        grid.swapValidMasks();
        grid.swapVelocityBuffers();
    }
}

int particles_num = 0;

void FluidSimulator::advance_flip_pic(FluidDomain &domain, float t_frame, float flip_pic_ratio = 0.98) {
    float t = 0;
    while (t < t_frame) {
        float timestep = compute_cfl(domain);
        if (timestep + t >= t_frame) { timestep = t_frame - t; }
        printf("Taking substep of size %f (to %0.3f%% of the frame)\n", timestep, 100 * (t + timestep) / t_frame);


        auto t_start = std::chrono::high_resolution_clock::now();
        domain.update(timestep);
        auto t_end = std::chrono::high_resolution_clock::now();
        auto update_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);

        t_start = std::chrono::high_resolution_clock::now();
        domain.constructFluidLevelSetFromMarkerParticles();
        t_end = std::chrono::high_resolution_clock::now();
        auto constructFluidLevelSetFromMarkerParticles_duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);

        t_start = std::chrono::high_resolution_clock::now();
        // reseeding(domain);
        t_end = std::chrono::high_resolution_clock::now();
        auto reseeding_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "before transfer");

        if (particles_num != domain.particleSet().size()) {
            particles_num = domain.particleSet().size();
            std::cout << "Number of particles: " << particles_num << std::endl;
        }

        t_start = std::chrono::high_resolution_clock::now();
        transfer_from_particles_to_grid(domain);
        t_end = std::chrono::high_resolution_clock::now();
        auto transfer_from_particles_to_grid_duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after transfer");

        domain.grid().validUFront()->tree().topologyUnion(domain.grid().velFront()->tree());
        domain.grid().validVFront()->tree().topologyUnion(domain.grid().velFront()->tree());
        domain.grid().validWFront()->tree().topologyUnion(domain.grid().velFront()->tree());
        t_start = std::chrono::high_resolution_clock::now();
        extrapolate_data(domain, 2);
        t_end = std::chrono::high_resolution_clock::now();
        auto extrapolate_data_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after extrapolation");

        domain.grid().updatePreviousVelocityBuffer();

        // Eulerian grid part START

        t_start = std::chrono::high_resolution_clock::now();
        add_forces(domain, timestep);
        t_end = std::chrono::high_resolution_clock::now();
        auto add_forces_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after forces");

        t_start = std::chrono::high_resolution_clock::now();
        project(domain, timestep);
        t_end = std::chrono::high_resolution_clock::now();
        auto project_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after project");
        t_start = std::chrono::high_resolution_clock::now();
        extrapolate_data(domain, 2);
        t_end = std::chrono::high_resolution_clock::now();
        extrapolate_data_duration += std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after extrapolation");

        t_start = std::chrono::high_resolution_clock::now();
        constrain_velocity(domain);
        t_end = std::chrono::high_resolution_clock::now();
        auto constrain_velocty_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
        print_velocity_field(domain.grid(), "after constrain");

        // Eulerian grid part END

        domain.grid().updateDiffBuffers();

        t_start = std::chrono::high_resolution_clock::now();
        transfer_from_grid_to_particles(domain, flip_pic_ratio);
        t_end = std::chrono::high_resolution_clock::now();
        auto transfer_from_grid_to_particles_duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);

        t_start = std::chrono::high_resolution_clock::now();
        domain.advectParticles(timestep);
        t_end = std::chrono::high_resolution_clock::now();
        auto advectParticles_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);


        printf(
            "update %.3f, construct %.3f, reseeding_duration %.3f, transfer grid %.3f, add forces %.3f, project %.3f, extrapolate %.3f, constrain_velocty_duration %.3f, transfer particles %.3f, advect particles %.3f\n",
            update_duration.count() / 1000.f, constructFluidLevelSetFromMarkerParticles_duration.count() / 1000.f,
            reseeding_duration.count() / 1000.f, transfer_from_particles_to_grid_duration.count() / 1000.f,
            add_forces_duration.count() / 1000.f, project_duration.count() / 1000.f,
            extrapolate_data_duration.count() / 1000.f, constrain_velocty_duration.count() / 1000.f,
            transfer_from_grid_to_particles_duration.count() / 1000.f, advectParticles_duration.count() / 1000.f);
        t += timestep;
    }
}

float FluidSimulator::compute_cfl(FluidDomain &domain) {
    auto &grid = domain.grid();
    openvdb::Vec3d minVector, maxVector;
    grid.velFront()->evalMinMax(minVector, maxVector);
    double max_vel = max(maxVector[0], maxVector[1], maxVector[2]);
    max_vel += sqrt(5 * max(domain.voxelSize(), domain.voxelSize()) * grav);
    return min(domain.voxelSize(), domain.voxelSize()) / (max_vel + epsilon);
}

void FluidSimulator::add_forces(FluidDomain &domain, float dt) {
    // gravity
    auto &grid = domain.grid();
    auto accessor = domain.fluidLevelSet().getAccessor();
    auto vel_accessor = grid.velFront()->getAccessor();
    for (auto iter = grid.velFront()->cbeginValueOn(); iter; ++iter) {
        auto coord = iter.getCoord();
        // Only add force to the liquid cells
        // if (accessor.getValue(coord) < 0) {
        if (VERTICAL_PLANE) {
            grid.setVelYHalfIndexed(vel_accessor, coord, grid.velHalfIndexed(vel_accessor, coord)[1] - grav * dt);
        }
        // }
    }
}

void FluidSimulator::project(FluidDomain &domain, float dt) {
    auto &grid = domain.grid();
    int system_size = 0;
    // Solver of linear system

    Eigen::setNbThreads(8);
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double, 1 >, Eigen::UpLoType::Lower | Eigen::UpLoType::Upper > solver;
    solver.setTolerance(1e-7);
    // Laplacian matrix
    Eigen::SparseMatrix<double, 1> alpha_matrix;

    openvdb::Int32Grid::Ptr fluid_indices = openvdb::createGrid<openvdb::Int32Grid>(-1);
    auto fluid_indices_accessor = fluid_indices->getAccessor();
    // Index all cells with fluid
    auto fluidAccessor = domain.fluidLevelSet().getAccessor();
    auto solidAccessor = domain.solidLevelSet().getAccessor();
    auto vel_front_accessor = grid.velFront()->getAccessor();
    auto vel_back_accessor = grid.velBack()->getAccessor();
    auto bbox = domain.fluidLevelSet().getLevelSet()->evalActiveVoxelBoundingBox();

    for (auto it = bbox.beginZYX(); it != bbox.endZYX(); ++it) {
        auto coord = (*it);
        if (fluidAccessor.getValue(coord) < 0) {
            fluid_indices_accessor.setValue(coord, system_size);
            system_size++;
        }
    }

    if (system_size == 0) { return; }
    alpha_matrix.resize(system_size, system_size);
    alpha_matrix.reserve(Eigen::VectorXi::Constant(system_size, 7));

    Eigen::VectorXd rhs(system_size);

    // Pressure solver
    // std::vector<double> rhs(system_size);
    // std::vector<double> pressure_grid(system_size);
    // SparseMatrixd alpha_matrix(system_size);
    // PCGSolver<double> solver;

    // alpha_matrix.zero();
    // rhs.assign(rhs.size(), 0);
    // pressure_grid.assign(pressure_grid.size(), 0);


    float scale_x = dt / (FluidDomain::density * sqr(domain.voxelSize()));
    float scale_y = dt / (FluidDomain::density * sqr(domain.voxelSize()));
    float scale_z = dt / (FluidDomain::density * sqr(domain.voxelSize()));

    // solid velocity zero for now
    float u_solid = 0;
    float v_solid = 0;
    float w_solid = 0;

    grid.uWeights()->clear();
    grid.vWeights()->clear();
    grid.wWeights()->clear();

    auto t_start = std::chrono::high_resolution_clock::now();
    // Resample solid level set to get values at bottom left corner of voxels
    auto resampledSolidLevelSet = openvdb::createLevelSet<openvdb::FloatGrid>(domain.voxelSize());
    openvdb::tools::GridTransformer transformer(
        domain.solidLevelSet().getLevelSet()->transformPtr()->baseMap()->getAffineMap()->getMat4()
        * resampledSolidLevelSet->transformPtr()->baseMap()->getAffineMap()->getMat4().inverse());

    transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(*domain.solidLevelSet().getLevelSet(),
                                                                              *resampledSolidLevelSet);
    auto border_mask = openvdb::tools::extractIsosurfaceMask(*domain.solidLevelSet().getLevelSet(), 0);
    openvdb::tools::dilateActiveValues(border_mask->tree(), 1, openvdb::tools::NN_FACE_EDGE_VERTEX);

    grid.uWeights()->tree().topologyUnion(border_mask->tree());
    grid.wWeights()->tree().topologyUnion(border_mask->tree());
    grid.vWeights()->tree().topologyUnion(border_mask->tree());

    openvdb::tree::IteratorRange<openvdb::BoolGrid::ValueOnCIter> range(border_mask->cbeginValueOn());
    tbb::parallel_for(range, [&](openvdb::tree::IteratorRange<openvdb::BoolGrid::ValueOnCIter> range) {
        auto targetAccessor = resampledSolidLevelSet->getAccessor();
        auto u_weights_accessor = grid.uWeights()->getAccessor();
        auto v_weights_accessor = grid.vWeights()->getAccessor();
        auto w_weights_accessor = grid.wWeights()->getAccessor();
        for (; range; ++range) {
            auto coord = range.iterator().getCoord();

            auto point_000 = targetAccessor.getValue(coord);
            auto point_001 = targetAccessor.getValue(coord.offsetBy(0, 0, 1));
            auto point_010 = targetAccessor.getValue(coord.offsetBy(0, 1, 0));
            auto point_011 = targetAccessor.getValue(coord.offsetBy(0, 1, 1));
            auto point_100 = targetAccessor.getValue(coord.offsetBy(1, 0, 0));
            auto point_101 = targetAccessor.getValue(coord.offsetBy(1, 0, 1));
            auto point_110 = targetAccessor.getValue(coord.offsetBy(1, 1, 0));

            u_weights_accessor.setValue(coord,
                                        1 - LevelSet::fraction_inside(point_000, point_010, point_001, point_011));
            v_weights_accessor.setValue(coord,
                                        1 - LevelSet::fraction_inside(point_000, point_001, point_100, point_101));
            w_weights_accessor.setValue(coord,
                                        1 - LevelSet::fraction_inside(point_000, point_010, point_100, point_110));
        }
    });
    auto t_end = std::chrono::high_resolution_clock::now();
    auto weights_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);


    t_start = std::chrono::high_resolution_clock::now();
    tbb::parallel_for(tbb::blocked_range<int>(_bbox.min()[0], _bbox.max()[0]), [&](tbb::blocked_range<int> i_range) {
        auto fluidAccessor = domain.fluidLevelSet().getAccessor();
        auto fluid_indices_accessor = fluid_indices->getAccessor();
        auto u_weights_accessor = grid.uWeights()->getAccessor();
        auto v_weights_accessor = grid.vWeights()->getAccessor();
        auto w_weights_accessor = grid.wWeights()->getAccessor();
        for (int i = i_range.begin(); i < i_range.end(); ++i) {
            for (int j = _bbox.min()[1]; j < _bbox.max()[1]; ++j) {
                for (int k = _bbox.min()[0]; k < _bbox.max()[0]; ++k) {
                    auto coord = openvdb::Coord(i, j, k);

                    if (fluid_indices_accessor.getValue(coord) != -1) {
                        int idx = fluid_indices_accessor.getValue(coord);
                        rhs[idx] = 0;
                        double diagonal = 0;
                        double theta = 0;
                        double term;

                        float phi_i_j_k = fluidAccessor.getValue(coord);
                        float phi_i_minus_1_j_k = fluidAccessor.getValue(coord.offsetBy(-1, 0, 0));
                        float phi_i_j_minus_1_k = fluidAccessor.getValue(coord.offsetBy(0, -1, 0));
                        float phi_i_j_k_minus_1 = fluidAccessor.getValue(coord.offsetBy(0, 0, -1));
                        float phi_i_plus_1_j_k = fluidAccessor.getValue(coord.offsetBy(1, 0, 0));
                        float phi_i_j_plus_1_k = fluidAccessor.getValue(coord.offsetBy(0, 1, 0));
                        float phi_i_j_k_plus_1 = fluidAccessor.getValue(coord.offsetBy(0, 0, 1));

                        term = u_weights_accessor.getValue(coord) * scale_x;
                        if (phi_i_minus_1_j_k < 0) {
                            alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(-1, 0, 0))) = -term;
                        }
                        theta = max(0.001, LevelSet::fraction_inside(phi_i_minus_1_j_k, phi_i_j_k));
                        diagonal += term * (1.f / theta);

                        rhs[idx] += u_weights_accessor.getValue(coord)
                                        * grid.velHalfIndexed(vel_front_accessor, coord)[0] / domain.voxelSize()
                                    + (1 - u_weights_accessor.getValue(coord)) * u_solid / domain.voxelSize();


                        term = u_weights_accessor.getValue(coord.offsetBy(1, 0, 0)) * scale_x;
                        if (phi_i_plus_1_j_k < 0) {
                            alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(1, 0, 0))) = -term;
                        }
                        theta = max(0.001, LevelSet::fraction_inside(phi_i_j_k, phi_i_plus_1_j_k));
                        diagonal += term * (1.f / theta);

                        rhs[idx] -=
                            u_weights_accessor.getValue(coord.offsetBy(1, 0, 0))
                                * grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(1, 0, 0))[0]
                                / domain.voxelSize()
                            + (1 - u_weights_accessor.getValue(coord.offsetBy(1, 0, 0))) * u_solid / domain.voxelSize();

                        term = v_weights_accessor.getValue(coord) * scale_y;
                        if (phi_i_j_minus_1_k < 0) {
                            alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(0, -1, 0))) = -term;
                        }
                        theta = max(0.001, LevelSet::fraction_inside(phi_i_j_minus_1_k, phi_i_j_k));
                        diagonal += term * (1.f / theta);

                        rhs[idx] += v_weights_accessor.getValue(coord)
                                        * grid.velHalfIndexed(vel_front_accessor, coord)[1] / domain.voxelSize()
                                    + (1 - v_weights_accessor.getValue(coord)) * v_solid / domain.voxelSize();

                        term = v_weights_accessor.getValue(coord.offsetBy(0, 1, 0)) * scale_y;
                        if (phi_i_j_plus_1_k < 0) {
                            alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(0, 1, 0))) = -term;
                        }
                        theta = max(0.001, LevelSet::fraction_inside(phi_i_j_k, phi_i_j_plus_1_k));
                        diagonal += term * (1.f / theta);

                        rhs[idx] -=
                            v_weights_accessor.getValue(coord.offsetBy(0, 1, 0))
                                * grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 1, 0))[1]
                                / domain.voxelSize()
                            + (1 - v_weights_accessor.getValue(coord.offsetBy(0, 1, 0))) * v_solid / domain.voxelSize();

                        term = w_weights_accessor.getValue(coord) * scale_z;
                        if (phi_i_j_k_minus_1 < 0) {
                            alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(0, 0, -1))) = -term;
                        }
                        theta = max(0.001, LevelSet::fraction_inside(phi_i_j_k_minus_1, phi_i_j_k));
                        diagonal += term * (1.f / theta);

                        rhs[idx] += w_weights_accessor.getValue(coord)
                                        * grid.velHalfIndexed(vel_front_accessor, coord)[2] / domain.voxelSize()
                                    + (1 - w_weights_accessor.getValue(coord)) * w_solid / domain.voxelSize();


                        term = w_weights_accessor.getValue(coord.offsetBy(0, 0, 1)) * scale_z;
                        if (phi_i_j_k_plus_1 < 0) {
                            alpha_matrix.insert(idx, fluid_indices_accessor.getValue(coord.offsetBy(0, 0, 1))) = -term;
                        }
                        theta = max(0.001, LevelSet::fraction_inside(phi_i_j_k, phi_i_j_k_plus_1));
                        diagonal += term * (1.f / theta);

                        rhs[idx] -=
                            w_weights_accessor.getValue(coord.offsetBy(0, 0, 1))
                                * grid.velHalfIndexed(vel_front_accessor, coord.offsetBy(0, 0, 1))[2]
                                / domain.voxelSize()
                            + (1 - w_weights_accessor.getValue(coord.offsetBy(0, 0, 1))) * w_solid / domain.voxelSize();
                        //

                        alpha_matrix.insert(idx, idx) = diagonal;
                    }
                    // }
                }
            }
        }
    });
    t_end = std::chrono::high_resolution_clock::now();
    auto system_build_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);

    t_start = std::chrono::high_resolution_clock::now();
    Eigen::VectorXd pressure_grid(system_size);
    solver.compute(alpha_matrix);
    pressure_grid = solver.solve(rhs);
    t_end = std::chrono::high_resolution_clock::now();
    auto solver_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
    // for (auto it = rhs.begin(); it != rhs.end(); ++it) { std::cout << *it << " "; }
    // std::cout << std::endl;
    // for (auto it = pressure_grid.begin(); it != pressure_grid.end(); ++it) { std::cout << *it << " "; }
    // std::cout << std::endl;

    t_start = std::chrono::high_resolution_clock::now();

    domain.grid().validUFront()->clear();
    domain.grid().validVFront()->clear();
    domain.grid().validWFront()->clear();

    auto u_valid_acc = domain.grid().validUFront()->getAccessor();
    auto v_valid_acc = domain.grid().validVFront()->getAccessor();
    auto w_valid_acc = domain.grid().validWFront()->getAccessor();

    grid.velBack()->clear();
    vel_back_accessor = grid.velBack()->getAccessor();
    auto u_weights_accessor = grid.uWeights()->getAccessor();
    auto v_weights_accessor = grid.vWeights()->getAccessor();
    auto w_weights_accessor = grid.wWeights()->getAccessor();
    for (auto it = _bbox.beginZYX(); it != _bbox.endZYX(); ++it) {
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
        auto newVel = grid.velHalfIndexed(vel_back_accessor, coord);

        if (u_weights_accessor.getValue(coord) > 0 && (phi_i_j_k < 0 || phi_i_minus_1_j_k < 0)) {
            float theta = 1;
            if (solidAccessor.getValue(coord) <= 0 || solidAccessor.getValue(coord.offsetBy(-1, 0, 0)) <= 0) {
                theta = max(0.001, LevelSet::fraction_inside(phi_i_j_k, phi_i_minus_1_j_k));
            }
            newVel[0] = prevVel[0] - dt / FluidDomain::density * (p - p_i_minus1) / domain.voxelSize() / theta;
            u_valid_acc.setValueOn(coord);
        }

        if (v_weights_accessor.getValue(coord) > 0 && (phi_i_j_k < 0 || phi_i_j_minus_1_k < 0)) {
            float theta = 1;
            if (solidAccessor.getValue(coord) <= 0 || solidAccessor.getValue(coord.offsetBy(0, -1, 0)) <= 0) {
                theta = max(0.001, LevelSet::fraction_inside(phi_i_j_k, phi_i_j_minus_1_k));
            }
            newVel[1] = prevVel[1] - dt / FluidDomain::density * (p - p_j_minus1) / domain.voxelSize() / theta;
            v_valid_acc.setValueOn(coord);
        }

        if (w_weights_accessor.getValue(coord) > 0 && (phi_i_j_k < 0 || phi_i_j_k_minus_1 < 0)) {
            float theta = 1;
            if (solidAccessor.getValue(coord) <= 0 || solidAccessor.getValue(coord.offsetBy(0, 0, -1)) <= 0) {
                theta = max(0.001, LevelSet::fraction_inside(phi_i_j_k, phi_i_j_k_minus_1));
            }
            newVel[2] = prevVel[2] - dt / FluidDomain::density * (p - p_k_minus1) / domain.voxelSize() / theta;
            w_valid_acc.setValueOn(coord);
        }
        grid.setVelHalfIndexed(vel_back_accessor, coord, newVel);
    }
    t_end = std::chrono::high_resolution_clock::now();
    auto update_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);

    printf("weights_duration %.3f, system_build_duration %.3f, solver_duration %.3f, update_duration %.3f \n",
           weights_duration.count() / 1000.f, system_build_duration.count() / 1000.f, solver_duration.count() / 1000.f,
           update_duration.count() / 1000.f);

    grid.swapVelocityBuffers();
}

void FluidSimulator::enforceDirichlet(FluidDomain &domain) {
    auto &grid = domain.grid();
    auto solidAccessor = domain.solidLevelSet().getAccessor();
    auto vel_accessor = grid.velFront()->getAccessor();
    for (auto it = _bbox.beginZYX(); it != _bbox.endZYX(); ++it) {
        auto coord = (*it);
        if ((solidAccessor.getValue(coord.offsetBy(-1, 0, 0)) <= 0 && grid.velHalfIndexed(vel_accessor, coord)[0] < 0)
            || (solidAccessor.getValue(coord) <= 0 && grid.velHalfIndexed(vel_accessor, coord)[0] > 0))
            grid.setVelXHalfIndexed(vel_accessor, coord, 0);

        if ((solidAccessor.getValue(coord.offsetBy(0, -1, 0)) <= 0 && grid.velHalfIndexed(vel_accessor, coord)[1] < 0)
            || (solidAccessor.getValue(coord) <= 0 && grid.velHalfIndexed(vel_accessor, coord)[1] > 0))
            grid.setVelYHalfIndexed(vel_accessor, coord, 0);

        if ((solidAccessor.getValue(coord.offsetBy(0, 0, -1)) <= 0 && grid.velHalfIndexed(vel_accessor, coord)[2] < 0)
            || (solidAccessor.getValue(coord) <= 0 && grid.velHalfIndexed(vel_accessor, coord)[2] > 0))
            grid.setVelZHalfIndexed(vel_accessor, coord, 0);
    }
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

    openvdb::tools::GridSampler<openvdb::Vec3DGrid::Accessor, openvdb::tools::StaggeredBoxSampler>
        vel_staggered_x_sampler(grid.velFront()->getAccessor(), grid.velFront()->transform());
    openvdb::tools::GridSampler<openvdb::Vec3DGrid::Accessor, openvdb::tools::StaggeredBoxSampler>
        vel_staggered_y_sampler(grid.velFront()->getAccessor(), grid.velFront()->transform());
    openvdb::tools::GridSampler<openvdb::Vec3DGrid::Accessor, openvdb::tools::StaggeredBoxSampler>
        vel_staggered_z_sampler(grid.velFront()->getAccessor(), grid.velFront()->transform());

    auto unionMask = openvdb::createGrid<openvdb::MaskGrid>();
    unionMask->tree().topologyUnion(grid.uWeights()->tree());
    unionMask->tree().topologyUnion(grid.vWeights()->tree());
    unionMask->tree().topologyUnion(grid.wWeights()->tree());

    for (auto iter = unionMask->beginValueOn(); iter; ++iter) {
        auto coord = iter.getCoord();
        openvdb::Vec3d newVel = vel_back_accessor.getValue(coord);
        if (u_weights_accessor.getValue(coord) == 0) {
            // auto normal = grad_staggered_x_accessor.getValue(coord).unit();
            auto normal = grad_staggered_x_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0.5, 0, 0)).unit();
            // auto vel = vel_staggered_x_accessor.getValue(coord);
            auto vel = vel_staggered_x_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0.5, 0, 0));
            auto perp_component = normal.dot(vel);
            auto projection = perp_component * normal;
            auto tangential = vel - projection;
            newVel[0] = tangential[0];
        }
        if (v_weights_accessor.getValue(coord) == 0) {
            // auto normal = grad_staggered_y_accessor.getValue(coord).unit();
            auto normal = grad_staggered_y_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0.5, 0)).unit();
            // auto vel = vel_staggered_y_accessor.getValue(coord);
            auto vel = vel_staggered_y_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0.5, 0));
            auto perp_component = normal.dot(vel);
            auto projection = perp_component * normal;
            auto tangential = vel - projection;
            newVel[1] = tangential[1];
        }
        if (w_weights_accessor.getValue(coord) == 0) {
            // auto normal = grad_staggered_z_accessor.getValue(coord).unit();
            auto normal = grad_staggered_z_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0, 0.5)).unit();
            // auto vel = vel_staggered_z_accessor.getValue(coord);
            auto vel = vel_staggered_z_sampler.isSample(coord.asVec3d() - openvdb::Vec3d(0, 0, 0.5));
            auto perp_component = normal.dot(vel);
            auto projection = perp_component * normal;
            auto tangential = vel - projection;
            newVel[2] = tangential[2];
        }
        grid.setVelHalfIndexed(vel_back_accessor, coord, newVel);
    }

    grid.swapVelocityBuffers();
}

double FluidSimulator::calculate_kernel_function(double x, double y, double z) {
    double val1 = calculate_quad_bspline(x);
    double val2 = calculate_quad_bspline(y);
    double val3 = calculate_quad_bspline(z);
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