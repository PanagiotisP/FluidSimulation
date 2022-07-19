#include "ParticleSet.h"

#include "util.h"

#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/PointAdvect.h>
#include <openvdb/tools/TopologyToLevelSet.h>
#include <tbb/parallel_for.h>

Particle::Particle() {}
Particle::Particle(openvdb::Vec3d pos, openvdb::Vec3d vel, float radius, float mass): _pos(pos), _vel(vel), _radius(radius), _mass(mass) {};
Particle::~Particle() {}

ParticleSet::ParticleSet() {};
ParticleSet::ParticleSet(openvdb::math::Transform::Ptr i2w_transform): _i2w_transform(i2w_transform) {};
ParticleSet::~ParticleSet() {};

void ParticleSet::advectAndEnsureOutsideObstacles(openvdb::Vec3dGrid::Ptr vel_grid, openvdb::Vec3fGrid::Ptr cpt_grid,
                                                  LevelSet &solid_level_set, float dt) {  
    // Probably costly way of extracting positions
    // Extract positions (needed by PointAdvect.advect())
    std::vector<openvdb::Vec3d> positions(_particles.size());
    tbb::parallel_for(tbb::blocked_range<int>(0, _particles.size()), [&](tbb::blocked_range<int> r) {
        for (int i = r.begin(); i < r.end(); ++i) { positions[i] = _i2w_transform->indexToWorld(_particles[i].pos()); }
    });
    
    // Use RK-3 method to advect positions
    openvdb::tools::PointAdvect<openvdb::Vec3dGrid, std::vector<openvdb::Vec3d>, 1> p_advect(*vel_grid);
    p_advect.setIntegrationOrder(3);
    p_advect.advect(positions, dt);
    
    tbb::parallel_for(tbb::blocked_range<int>(0, _particles.size()), [&](tbb::blocked_range<int> r) {
        std::shared_ptr<LevelSet::Sampler> solidSampler(solid_level_set.getSampler(solid_level_set.getLevelSet()->getAccessor()));
        // LevelSet::Sampler solidSampler(
            // solid_level_set->getAccessor(), solid_level_set->transform());
        openvdb::tools::GridSampler<openvdb::Vec3fGrid::Accessor, openvdb::tools::BoxSampler> cptSampler(
            cpt_grid->getAccessor(), cpt_grid->transform());
        for (int i = r.begin(); i < r.end(); ++i) {
            if (solidSampler->wsSample(positions[i]) <= 0)
                _particles[i].setPosition(_i2w_transform->worldToIndex(cptSampler.wsSample(positions[i])));
            else
                _particles[i].setPosition(_i2w_transform->worldToIndex(positions[i]));
        }
    });
}