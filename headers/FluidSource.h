#pragma once

#include "LevelSet.h"
#include "MacGrid.h"

#include <openvdb/openvdb.h>

class FluidSource {
public:
    FluidSource();
    FluidSource(LevelSet spawning_region, openvdb::Vec3d vel, int particle_generation_rate, LevelSet &solidLevelSet, float max_spawning_duration);
    ~FluidSource();

    // Updates this fluid source by spawning new particles.
    // Reduces spawning duration left
    void update(MacGrid &grid, ParticleSet &particle_set, float dt);

    // Resets the fluid source, making it spawn again for the full duration
    void reset() {
        _spawning_duration_left = _max_spawning_duration;
    }

    inline LevelSet &spawningRegion() { return _spawning_region; };

private:
    openvdb::Vec3d _vel;

    LevelSet _spawning_region;
    float _spawning_duration_left;
    float _max_spawning_duration;
    float _particle_generation_rate;
};
