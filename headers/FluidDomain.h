#pragma once
#include "LevelSet.h"
#include "MacGrid.h"
#include "ParticleSet.h"


#define ZERO_CELCIUS 273

class FluidSource {
public:
    FluidSource(LevelSet spawning_region, float vel_x, float vel_y, float dx, float dy, int particle_generation_rate,
                float rate_t, float rate_c);
    FluidSource(int size_x, int size_y, float length_x, float length_y);
    ~FluidSource();
    void update(MacGrid &grid, ParticleSet &particle_set, float dt);

    inline LevelSet &spawningRegion() { return spawning_region; };

private:
    float dx;
    float dy;

    float _vel_x;
    float _vel_y;

    float rate_t;
    float rate_c;


    LevelSet spawning_region;
    float particle_generation_rate;
};

class FluidDomain {
public:
    FluidDomain(int size_x, int size_y, float length_x, float length_y);
    ~FluidDomain();

    // static for now
    static constexpr float density = 1;
    static constexpr float rate_t = 10;
    static constexpr float rate_c = 10;
    static constexpr float diffuse_cosnt_t = 1;
    static constexpr float ambient_temp = ZERO_CELCIUS + 21;
    static constexpr float ambient_concentration = 0.05;
    static constexpr float diffuse_rate_temperature = 0.1;
    static constexpr float diffuse_rate_concentration = 0;

    void addFluidSource(FluidSource fluid_source);
    void removeFluidSource(int i);
    void update(float dt);
    void clearParticleSet();
    void clearFluidSources();

    void construct_level_set_from_marker_particles(LevelSet &level_set);

    inline MacGrid &grid() { return _grid; };
    inline LevelSet &fluidLevelSet() { return fluid_level_set; };
    inline LevelSet &solidLevelSet() { return solid_level_set; };
    inline void setFluidLevelSet(LevelSet &l) { fluid_level_set = l; };
    inline void setSolidLevelSet(LevelSet &l) { solid_level_set = l; };
    inline ParticleSet &particleSet() { return particle_set; };
    inline std::vector<FluidSource> &fluidSources() { return fluid_sources; };

    void advectParticles(float dt);

    void classifyCells();

    float radius;
private:
    MacGrid _grid;
    LevelSet fluid_level_set;
    LevelSet solid_level_set;
    ParticleSet particle_set;
    std::vector<FluidSource> fluid_sources;
};