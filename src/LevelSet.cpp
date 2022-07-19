#include "LevelSet.h"

#include "util.h"

#include <Eigen/Geometry>
#include <functional>
#include <iomanip>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetAdvect.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/tools/VelocityFields.h>

LevelSet::~LevelSet() {}
LevelSet::LevelSet(openvdb::FloatGrid::Ptr level_set):
 _level_set(level_set), _accessor(level_set->tree()), _sampler(_accessor, _level_set->transform()) { }

void LevelSet::unionOfBalls(ParticleSet &points, float voxel_size) {
    _level_set->clear();
    openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*_level_set);
    raster.setRmin(points.getMinRadius());
    raster.rasterizeSpheres(points);
    raster.finalize(true);
}

float LevelSet::valueInterpolatedI(LevelSet::Sampler &sampler, float x, float y, float z) {
    return sampler.isSample(openvdb::Vec3R(x, y, z));
}
float LevelSet::valueInterpolatedI(std::shared_ptr<LevelSet::Sampler> &sampler, openvdb::Vec3R isPoint) {
    return sampler->isSample(isPoint);
}

float LevelSet::valueInterpolatedW(LevelSet::Sampler &sampler, float x, float y, float z) {
    return sampler.wsSample(openvdb::Vec3R(x, y, z));
}
float LevelSet::valueInterpolatedW(LevelSet::Sampler &sampler, openvdb::Vec3R wsPoint) {
    return sampler.wsSample(wsPoint);
}

void LevelSet::unionLevelSetInPlace(LevelSet &level_set) {
    _level_set = openvdb::tools::csgUnionCopy(*_level_set, *level_set.getLevelSet());
}
void LevelSet::intersectionLevelSetInPlace(LevelSet &level_set) {
    _level_set = openvdb::tools::csgIntersectionCopy(*_level_set, *level_set.getLevelSet());
}
void LevelSet::differenceLevelSetInPlace(LevelSet &level_set) {
    _level_set = openvdb::tools::csgDifferenceCopy(*_level_set, *level_set.getLevelSet());
}

void LevelSet::advect(openvdb::Vec3dGrid::Ptr vel, float dt) {
    openvdb::tools::DiscreteField<openvdb::Vec3dGrid, openvdb::tools::StaggeredBoxSampler> vel_field(*vel);
    openvdb::tools::LevelSetAdvection<
        openvdb::FloatGrid, openvdb::tools::DiscreteField<openvdb::Vec3dGrid, openvdb::tools::StaggeredBoxSampler>>
        l(*_level_set, vel_field);
    std::cout << l.advect(0, dt) << std::endl;
}

void LevelSet::invert() {
    for (auto iter = _level_set->tree().beginRootChildren(); iter; ++iter) {
        iter->negate();
    }
}

LevelSet unionLevelSet(const LevelSet &level_set_a, const LevelSet &level_set_b) {
    auto level_set = openvdb::tools::csgUnionCopy(*level_set_a.getLevelSet(), *level_set_b.getLevelSet());
    return LevelSet(level_set);
}
LevelSet intersectionLevelSet(const LevelSet &level_set_a, const LevelSet &level_set_b) {
    auto level_set = openvdb::tools::csgIntersectionCopy(*level_set_a.getLevelSet(), *level_set_b.getLevelSet());
    return LevelSet(level_set);
}
LevelSet differenceLevelSet(const LevelSet &level_set_a, const LevelSet &level_set_b) {
    auto level_set = openvdb::tools::csgDifferenceCopy(*level_set_a.getLevelSet(), *level_set_b.getLevelSet());
    return LevelSet(level_set);
}
