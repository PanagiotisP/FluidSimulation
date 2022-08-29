#pragma once
#include "FluidSource.h"
#include "LevelSet.h"
#include "MacGrid.h"
#include "ParticleSet.h"
#include "SolidObject.h"

#include <openvdb/openvdb.h>
#define ZERO_CELCIUS 273

class FluidDomain {
public:
    FluidDomain(openvdb::FloatGrid::Ptr boundary,
                openvdb::math::Transform::Ptr i2w_transform, float voxel_size);
    ~FluidDomain();

    // static for now
    // Density in kg/m^3
    static const float density;

    void addFluidSource(std::shared_ptr<FluidSource> fluid_source);
    void removeFluidSource(int i);

    void addSolidObj(std::shared_ptr<SolidObject> solid_obj);
    void removeSolidObj(int i);
    std::shared_ptr<SolidObject> getSolidObj(int i);

    // Updates the domain by calling the update method of each fluid source
    void update(float dt);
    void clearParticleSet();
    void clearFluidSources();
    void clearSolidObjs();

    void constructFluidLevelSetFromMarkerParticles();

    inline MacGrid &grid() { return _grid; };
    inline LevelSet &fluidLevelSet() { return _fluid_level_set; };
    inline LevelSet &solidLevelSet() { return _solid_level_set; };
    inline std::vector<std::shared_ptr<SolidObject>> &solidObjs() { return _solid_objs; }
    // inline void setFluidLevelSet(LevelSet &l) { fluid_level_set = l; };
    // inline void setSolidLevelSet(LevelSet &l) { solid_level_set = l; };
    inline ParticleSet &particleSet() { return particle_set; };
    inline std::vector<std::shared_ptr<FluidSource>> &fluidSources() { return _fluid_sources; };
    inline float voxelSize() { return _voxel_size; };

    void advectParticles(float dt);

    openvdb::math::Transform::Ptr _i2w_transform;

private:
    MacGrid _grid;
    LevelSet _fluid_level_set;
    LevelSet _solid_level_set;
    ParticleSet particle_set;
    std::vector<std::shared_ptr<FluidSource>> _fluid_sources;
    std::vector<std::shared_ptr<SolidObject>> _solid_objs;
    float _voxel_size;
};