#include "FluidDomain.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/FastSweeping.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <random>

FluidDomain::FluidDomain(openvdb::FloatGrid::Ptr boundary,
                         openvdb::math::Transform::Ptr i2w_transform, float voxel_size):
 _grid(i2w_transform, voxel_size),
 _solid_level_set(boundary),
 _fluid_level_set(openvdb::createLevelSet<openvdb::FloatGrid>(voxel_size)),
 particle_set(i2w_transform), _i2w_transform(i2w_transform), _voxel_size(voxel_size) {
    _fluid_level_set.getLevelSet()->setTransform(_i2w_transform);
    _solid_level_set.getLevelSet()->setTransform(_i2w_transform);
 }
FluidDomain::~FluidDomain() {}

void FluidDomain::update(float dt) {
    for (auto fluid_source : _fluid_sources) {
        fluid_source->update(_grid, particle_set, dt);
    }
}

void FluidDomain::addFluidSource(std::shared_ptr<FluidSource> fluid_source) { _fluid_sources.push_back(fluid_source); }

void FluidDomain::removeFluidSource(int i) {
    if (_fluid_sources.size() >= i + 1) _fluid_sources.erase(_fluid_sources.begin() + i);
}

void FluidDomain::addSolidObj(std::shared_ptr<SolidObject> solid_obj) { _solid_objs.push_back(solid_obj); }

void FluidDomain::removeSolidObj(int i) {
    if (_solid_objs.size() >= i + 1) _solid_objs.erase(_solid_objs.begin() + i);
}
std::shared_ptr<SolidObject> FluidDomain::getSolidObj(int i) {
    return _solid_objs[i];
}

void FluidDomain::clearParticleSet() { particle_set.clear(); }
void FluidDomain::clearFluidSources() { _fluid_sources.clear(); }
void FluidDomain::clearSolidObjs() { _solid_objs.clear(); }

void FluidDomain::constructFluidLevelSetFromMarkerParticles() {
    if (particle_set.size() > 0) {
        _fluid_level_set.unionOfBalls(particle_set, _voxel_size);
        auto fluidGridTopology = openvdb::createGrid<openvdb::MaskGrid>();
        fluidGridTopology->tree().topologyUnion(_fluid_level_set.getLevelSet()->tree());
        _fluid_level_set.getLevelSet()->tree().combineExtended(_solid_level_set.getLevelSet()->deepCopy()->tree(),
            // *openvdb::tools::csgUnionCopy(_solid.levelSet().getLevelSet()->deepCopy()->tree(),
                                        //  solid_level_set.getLevelSet()->deepCopy()->tree()),
            IntersectSolidFluidLevelSets::intersect, true);
        _fluid_level_set.getLevelSet()->tree().topologyIntersection(fluidGridTopology->tree());
    }
}

void FluidDomain::advectParticles(float dt) {
    LevelSet solid_objs_level_set(openvdb::createLevelSet<openvdb::FloatGrid>(_voxel_size));
    solid_objs_level_set.getLevelSet()->setTransform(_i2w_transform);

    for(auto solid_obj : _solid_objs) {
        solid_objs_level_set.unionLevelSetInPlace(solid_obj->levelSet());
    }

    auto combined_solid_level_set = unionLevelSet(_solid_level_set, solid_objs_level_set);
    auto cpt_grid = openvdb::tools::cpt(*combined_solid_level_set.getLevelSet());
    openvdb::tools::VelocityIntegrator<openvdb::Vec3dGrid, 1> rk_integrator(*_grid.velFront());
    particle_set.advectAndEnsureOutsideObstacles(_grid.velFront(), cpt_grid, combined_solid_level_set, dt);
}
