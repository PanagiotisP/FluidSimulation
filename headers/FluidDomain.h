#pragma once
#include "LevelSet.h"
#include "MacGrid.h"
#include "ParticleSet.h"

#include <openvdb/openvdb.h>
#define ZERO_CELCIUS 273

class FluidSource {
public:
    FluidSource();
    FluidSource(LevelSet spawning_region, openvdb::Vec3d vel, int particle_generation_rate);
    ~FluidSource();
    void update(MacGrid &grid, ParticleSet &particle_set, float dt);

    inline LevelSet &spawningRegion() { return spawning_region; };

private:
    openvdb::Vec3d _vel;

    LevelSet spawning_region;
    float particle_generation_rate;
};

class FluidDomain {
public:
    FluidDomain(openvdb::FloatGrid::Ptr solid_level_set,
                openvdb::FloatGrid::Ptr fluid_level_set, openvdb::math::Transform::Ptr i2w_transform, float voxel_size);
    ~FluidDomain();

    // static for now
    static constexpr float density = 1;

    void addFluidSource(FluidSource fluid_source);
    void removeFluidSource(int i);
    void update(float dt);
    void clearParticleSet();
    void clearFluidSources();

    void constructFluidLevelSetFromMarkerParticles();

    inline MacGrid &grid() { return _grid; };
    inline LevelSet &fluidLevelSet() { return fluid_level_set; };
    inline LevelSet &solidLevelSet() { return solid_level_set; };
    // inline void setFluidLevelSet(LevelSet &l) { fluid_level_set = l; };
    // inline void setSolidLevelSet(LevelSet &l) { solid_level_set = l; };
    inline ParticleSet &particleSet() { return particle_set; };
    inline std::vector<FluidSource> &fluidSources() { return fluid_sources; };
    inline float voxelSize() { return voxel_size; };

    void advectParticles(float dt);

    void classifyCells();

    openvdb::math::Transform::Ptr i2w_transform;

    struct IntersectSolidFluidLevelSets {
        static inline void intersect(openvdb::CombineArgs<float>& args) {
            // Args a fluid, Args b solid
            // If inside the solid then either the fluid is not touching (max = fluid_val), or it is touching and is -solid distance away
            float result;
            if (args.b() <= 0)
                result = max(args.a(), -args.b());
            else
                result = args.a();
            args.setResult(result);
            args.setResultIsActive(args.aIsActive());
        }
    };

private:
    MacGrid _grid;
    LevelSet fluid_level_set;
    LevelSet solid_level_set;
    ParticleSet particle_set;
    std::vector<FluidSource> fluid_sources;
    float voxel_size;
};