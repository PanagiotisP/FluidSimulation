#include "ParticleSet.h"

#include "util.h"

#include <openvdb/tools/Interpolation.h>
#include <tbb/parallel_for.h>

Particle::Particle() {}
Particle::Particle(openvdb::Vec3d pos, openvdb::Vec3d vel): _pos(pos), _vel(vel) {};
Particle::~Particle() {}

ParticleSet::ParticleSet() {};
ParticleSet::ParticleSet(openvdb::math::Transform::Ptr i2w_transform): i2w_transform(i2w_transform) {};
ParticleSet::~ParticleSet() {};

void ParticleSet::addParticle(std::unique_ptr<Particle> &p) { particles.push_back(std::move(p)); }
void ParticleSet::removeParticle(ParticleSet::iterator p) { particles.erase(p); }
void ParticleSet::advect(float dt) {
    for (auto it = particles.begin(); it != particles.end(); ++it) (*it)->advect(dt, i2w_transform);
}

void ParticleSet::advectAndEnsureOutsideObstacles(LevelSet &solid_level_set, openvdb::Vec3fGrid::Ptr cpt_grid,
                                                  float dt) {
    // LevelSet::BoxSampler sampler(solid_level_set.getAccessor(), solid_level_set.getLevelSet()->transform());
    // openvdb::tools::GridSampler<openvdb::Vec3fGrid::Accessor, openvdb::tools::BoxSampler> cpt_sampler(
    //     cpt_grid->getAccessor(), cpt_grid->transform());

    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> sampler(*solid_level_set.getLevelSet());
    openvdb::tools::GridSampler<openvdb::Vec3fGrid, openvdb::tools::BoxSampler> cpt_sampler(*cpt_grid);

    tbb::parallel_for(tbb::blocked_range<int>(0, particles.size()), [&](tbb::blocked_range<int> &range) {
        for (int i = range.begin(); i < range.end(); ++i) {
            particles[i]->advect(dt, i2w_transform);
            if (sampler.isSample(particles[i]->pos()) < 0) {
                particles[i]->setPosition(
                    i2w_transform->worldToIndex(openvdb::Vec3d(cpt_sampler.isSample(particles[i]->pos()))));
            }
        }
    });
}