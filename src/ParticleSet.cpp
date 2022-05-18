#include "ParticleSet.h"

#include "util.h"

#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/PointAdvect.h>
#include <tbb/parallel_for.h>

Particle::Particle() {}
Particle::Particle(openvdb::Vec3d pos, openvdb::Vec3d vel): _pos(pos), _vel(vel) {};
Particle::~Particle() {}

ParticleSet::ParticleSet() {};
ParticleSet::ParticleSet(openvdb::math::Transform::Ptr i2w_transform): i2w_transform(i2w_transform) {};
ParticleSet::~ParticleSet() {};

void ParticleSet::addParticle(const Particle &p) { particles.push_back(p); }
void ParticleSet::removeParticle(ParticleSet::iterator p) { particles.erase(p); }
void ParticleSet::advectAndEnsureOutsideObstacles(openvdb::Vec3dGrid::Ptr vel_grid, openvdb::Vec3fGrid::Ptr cpt_grid,
                                                  openvdb::FloatGrid::Ptr solid_level_set, float dt) {  
    // Probably costly way of extracting positions
    std::vector<openvdb::Vec3d> positions(particles.size());
    tbb::parallel_for(tbb::blocked_range<int>(0, particles.size()), [&](tbb::blocked_range<int> r) {
        for (int i = r.begin(); i < r.end(); ++i) { positions[i] = i2w_transform->indexToWorld(particles[i].pos()); }
    });
    openvdb::tools::PointAdvect<openvdb::Vec3dGrid, std::vector<openvdb::Vec3d>, 1> p_advect(*vel_grid);
    p_advect.setIntegrationOrder(3);
    p_advect.advect(positions, dt);
    tbb::parallel_for(tbb::blocked_range<int>(0, particles.size()), [&](tbb::blocked_range<int> r) {
        LevelSet::Sampler solidSampler(
            solid_level_set->getAccessor(), solid_level_set->transform());
        openvdb::tools::GridSampler<openvdb::Vec3fGrid::Accessor, openvdb::tools::BoxSampler> cptSampler(
            cpt_grid->getAccessor(), cpt_grid->transform());
        for (int i = r.begin(); i < r.end(); ++i) {
            if (solidSampler.wsSample(positions[i]) <= 0)
                particles[i].setPosition(i2w_transform->worldToIndex(cptSampler.wsSample(positions[i])));
            else
                particles[i].setPosition(i2w_transform->worldToIndex(positions[i]));
        }
    });
}