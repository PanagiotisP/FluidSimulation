#include "FluidDomain.h"

FluidSource::FluidSource() {};
FluidSource::FluidSource(float vel_x, float vel_y, float dx, float dy, int particle_generation_rate, BBox<float> area,
                         float rate_t = FluidDomain::rate_t, float rate_c = FluidDomain::rate_c):
 _vel_x(vel_x),
 _vel_y(vel_y), dx(dx), dy(dy), particle_generation_rate(particle_generation_rate), spawning_region(area),
 rate_t(rate_t), rate_c(rate_c) {}

FluidSource::~FluidSource() {}

void FluidSource::update(ParticleSet &particle_set, float dt) {
    float rate = dt * particle_generation_rate;
    int spawn_counter = static_cast<int>(rate);
    float probability_threshold = rate - static_cast<float>(spawn_counter);
    float particle_temperature = dt * rate_t;
    float particle_concentration = dt * rate_c;
    int i_start = static_cast<int>(spawning_region.x_min / dx);
    int j_start = static_cast<int>(spawning_region.y_min / dy);
    int i_end = static_cast<int>(spawning_region.x_max / dx);
    int j_end = static_cast<int>(spawning_region.y_max / dy);

    for (int j = j_start; j <= j_end; ++j) {
        for (int i = i_start; i <= i_end; ++i) {
            int grid_particles = spawn_counter;
            if (std::rand() / static_cast<float>(RAND_MAX) < probability_threshold) grid_particles++;
            while (grid_particles-- > 0) {
                auto p = std::make_unique<Particle>((i + std::rand() / static_cast<float>(RAND_MAX)) * dx,
                                                    (j + std::rand() / static_cast<float>(RAND_MAX)) * dy, _vel_x,
                                                    _vel_y, particle_temperature, particle_concentration);
                particle_set.addParticle(p);
            }
        }
    }
}


FluidDomain::FluidDomain(int size_x, int size_y, float length_x, float length_y):
 _grid(size_x, size_y, length_x, length_y, FluidDomain::ambient_temp, 0),
 level_set(size_x, size_y, length_x, length_y) {}
FluidDomain::~FluidDomain() {}

void FluidDomain::update(float dt) {
    for (auto it = fluid_sources.begin(); it != fluid_sources.end(); ++it) { it->update(particle_set, dt); }
}
void FluidDomain::removeFluidSource(int i) { fluid_sources.erase(fluid_sources.begin() + i); }

void FluidDomain::clearParticleSet() { particle_set.clear(); }
void FluidDomain::clearFluidSources() { fluid_sources.clear(); }

void FluidDomain::advectParticles(float dt) { particle_set.advectAndEnsureOutsideObstacles(_grid, dt); }

void FluidDomain::addFluidSource(FluidSource fluid_source) { fluid_sources.push_back(fluid_source); }

void FluidDomain::classifyCells() {}
void FluidDomain::classifyCells(ParticleSet &particle_set) {}