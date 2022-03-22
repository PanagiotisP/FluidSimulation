#pragma once
#include "LevelSet.h"
#include "MacGrid.h"
#include "ParticleSet.h"


#define ZERO_CELCIUS 273

class FluidSource {
public:
    FluidSource();
    FluidSource(float vel_x, float vel_y,float dx, float dy, int particle_generation_rate, BBox<float> area, float rate_t, float rate_c);
    ~FluidSource();
    void update(ParticleSet &particle_set, float dt);

private:
    float dx;
    float dy; 

    float _vel_x;
    float _vel_y;

    float rate_t;
    float rate_c;


    BBox<float> spawning_region;
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

    inline MacGrid &grid() {return _grid;};
    inline LevelSet &levelSet(){return level_set;};
    inline ParticleSet &particleSet(){return particle_set;};
    inline std::vector<FluidSource> &fluidSources() { return fluid_sources; };

    void advectParticles(float dt);

    void classifyCells();
    void classifyCells(ParticleSet &particle_set);

private:
    MacGrid _grid;
    LevelSet level_set;
    ParticleSet particle_set;
    std::vector<FluidSource> fluid_sources;
};