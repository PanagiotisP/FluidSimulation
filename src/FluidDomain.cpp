#include "FluidDomain.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/FastSweeping.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <random>

FluidSource::FluidSource(): spawning_region(openvdb::createLevelSet<openvdb::FloatGrid>()) {}
FluidSource::FluidSource(LevelSet spawning_region, openvdb::Vec3d vel, int particle_generation_rate):
 _vel(vel), particle_generation_rate(particle_generation_rate), spawning_region(spawning_region) {}

FluidSource::~FluidSource() {}

void FluidSource::update(MacGrid &grid, ParticleSet &particle_set, float dt) {
    float rate = dt * particle_generation_rate;
    int spawn_counter = static_cast<int>(rate);
    float probability_threshold = rate - static_cast<float>(spawn_counter);
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(-0.5, 0.5);

    auto coordBBox = spawning_region.getActiveCoordBBox();
    auto accessor = spawning_region.getAccessor();
    LevelSet::BoxSampler sampler(accessor, spawning_region.getLevelSet()->transform());

    MacGrid::BoxSampler vel_sampler(grid.velFront()->getAccessor(), spawning_region.getLevelSet()->transform());

    assert(coordBBox.min() < coordBBox.max());
    for (auto it = coordBBox.beginZYX(); it != coordBBox.endZYX(); ++it) {
        auto cellCenterCoord = (*it);
        if (accessor.getValue(cellCenterCoord) < 0) {
            int grid_particles = spawn_counter;
            if (dis(gen) < probability_threshold) grid_particles++;
            while (grid_particles > 0) {
                openvdb::Vec3d new_pos = cellCenterCoord.asVec3d() + openvdb::Vec3d(dis(gen), dis(gen), dis(gen));
                if (spawning_region.valueInterpolatedI(sampler, new_pos) < 0) {
                    openvdb::Vec3d new_vel = grid.velInterpolatedI(vel_sampler, new_pos) + _vel;
                    Particle p(new_pos, new_vel);
                    particle_set.addParticle(p);
                    grid_particles--;
                }
            }
        }
    }
}

FluidDomain::FluidDomain(openvdb::FloatGrid::Ptr solid_level_set, openvdb::FloatGrid::Ptr fluid_level_set,
                         openvdb::math::Transform::Ptr i2w_transform, float voxel_size):
 _grid(i2w_transform),
 solid_level_set(solid_level_set), fluid_level_set(fluid_level_set), particle_set(i2w_transform),
 i2w_transform(i2w_transform), voxel_size(voxel_size) {
    particle_set.setRadius(1.01f * voxel_size * sqrt(3) / 2.f);
}
FluidDomain::~FluidDomain() {}

void FluidDomain::update(float dt) {
    for (auto it = fluid_sources.begin(); it != fluid_sources.end(); ++it) { it->update(_grid, particle_set, dt); }
}

void FluidDomain::removeFluidSource(int i) {
    if (fluid_sources.size() >= i + 1) fluid_sources.erase(fluid_sources.begin() + i);
}

void FluidDomain::clearParticleSet() { particle_set.clear(); }
void FluidDomain::clearFluidSources() { fluid_sources.clear(); }

void FluidDomain::constructFluidLevelSetFromMarkerParticles() {
    if (particle_set.size() > 0) {
        fluid_level_set.construct_from_points(particle_set, voxel_size);
        auto fluidGridTopology = openvdb::createGrid<openvdb::MaskGrid>();
        fluidGridTopology->tree().topologyUnion(fluid_level_set.getLevelSet()->tree());
        fluid_level_set.getLevelSet()->tree().combineExtended(solid_level_set.getLevelSet()->deepCopy()->tree(),
                                                              IntersectSolidFluidLevelSets::intersect, true);
        fluid_level_set.getLevelSet()->tree().topologyIntersection(fluidGridTopology->tree());
    }
}

void FluidDomain::advectParticles(float dt) {
    auto cpt_grid = openvdb::tools::cpt(*solid_level_set.getLevelSet());
    openvdb::tools::VelocityIntegrator<openvdb::Vec3dGrid, 1> rk_integrator(*_grid.velFront());
    particle_set.advectAndEnsureOutsideObstacles(_grid.velFront(), cpt_grid, solid_level_set.getLevelSet(), dt);
}

void FluidDomain::addFluidSource(FluidSource fluid_source) { fluid_sources.push_back(fluid_source); }
