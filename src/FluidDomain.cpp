#include "FluidDomain.h"

#include <random>

FluidSource::FluidSource(int size_x, int size_y, float length_x, float length_y):
 spawning_region(size_x, size_y, length_x, length_y) {};
FluidSource::FluidSource(LevelSet spawning_region, float vel_x, float vel_y, float dx, float dy,
                         int particle_generation_rate, float rate_t = FluidDomain::rate_t,
                         float rate_c = FluidDomain::rate_c):
 _vel_x(vel_x),
 _vel_y(vel_y), dx(dx), dy(dy), particle_generation_rate(particle_generation_rate), rate_t(rate_t), rate_c(rate_c),
 spawning_region(spawning_region) {}

FluidSource::~FluidSource() {}

void FluidSource::update(MacGrid &grid, ParticleSet &particle_set, float dt) {
    float rate = dt * particle_generation_rate;
    int spawn_counter = static_cast<int>(rate);
    float probability_threshold = rate - static_cast<float>(spawn_counter);
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0, 1.0);
    for (int j = 0; j < spawning_region.sizeY(); ++j) {
        for (int i = 0; i < spawning_region.sizeX(); ++i) {
            if (spawning_region(i, j) < 0) {
                int grid_particles = spawn_counter;
                if (dis(gen) < probability_threshold) grid_particles++;
                while (grid_particles > 0) {
                    float new_x = (i + dis(gen)) * dx;
                    float new_y = (j + dis(gen)) * dy;
                    if (spawning_region.valueInterpolated(new_x, new_y) < 0) {
                        float current_temperature = grid.temperatureInterpolated(new_x, new_y);
                        float current_concentration = grid.concentrationInterpolated(new_x, new_y);

                        float new_vel_x = grid.velXInterpolated(new_x, new_y) + _vel_x;
                        float new_vel_y = grid.velYInterpolated(new_x, new_y) + _vel_y;

                        float new_temp = grid.temperatureInterpolated(new_x, new_y) + rate_t * dt;
                        float new_conc = clamp(grid.concentrationInterpolated(new_x, new_y) + rate_c * dt, 0.f, 1.f);
                        auto p = std::make_unique<Particle>(new_x, new_y, new_vel_x, new_vel_y, new_temp, new_conc);
                        particle_set.addParticle(p);
                        grid_particles--;
                    }
                }
            }
        }
    }
}

FluidDomain::FluidDomain(int size_x, int size_y, float length_x, float length_y):
 _grid(size_x, size_y, length_x, length_y, FluidDomain::ambient_temp, 0),
 fluid_level_set(size_x, size_y, length_x, length_y), solid_level_set(size_x, size_y, length_x, length_y) {
    //  little less than grid size
    radius = 1.01f * sqrt(sqr(_grid.deltaX()) + sqr(_grid.deltaY())) / 2.f;
}
FluidDomain::~FluidDomain() {}

void FluidDomain::update(float dt) {
    for (auto it = fluid_sources.begin(); it != fluid_sources.end(); ++it) {
        it->update(_grid, particle_set, dt);
        fluid_level_set.union_level_set(it->spawningRegion());
    }
}

void FluidDomain::removeFluidSource(int i) {
    if (fluid_sources.size() >= i + 1) fluid_sources.erase(fluid_sources.begin() + i);
}

void FluidDomain::clearParticleSet() { particle_set.clear(); }
void FluidDomain::clearFluidSources() { fluid_sources.clear(); }

void FluidDomain::construct_level_set_from_marker_particles(LevelSet& level_set) {
    std::vector<Vec2f> marker_points;
    marker_points.reserve(particle_set.size());
    for (auto it = particle_set.begin(); it != particle_set.end(); ++it) {
        auto &particle = *it;
        marker_points.push_back(Vec2f(particle->posX(), particle->posY()));
    }
    level_set.construct_from_points(marker_points);
    level_set.add_to_set(-radius);
    LevelSet inverted(solid_level_set.invert());
    level_set.intersect_level_set(inverted);
    level_set.redistance();
    
}

void FluidDomain::advectParticles(float dt) { particle_set.advectAndEnsureOutsideObstacles(solid_level_set, _grid, dt); }

void FluidDomain::addFluidSource(FluidSource fluid_source) { fluid_sources.push_back(fluid_source); }

void FluidDomain::classifyCells() {
    _grid.clearCellTypeBuffer();
    for (auto it = particle_set.begin(); it != particle_set.end(); ++it) {
        auto &p = *it;
        int i = clamp(static_cast<int>(floor(p->posX() / _grid.deltaX())), 0, _grid.sizeX() - 1);
        int j = clamp(static_cast<int>(floor(p->posY() / _grid.deltaY())), 0, _grid.sizeY() - 1);

        _grid.setCellType(i, j, LIQUID);
    }

    // Reset border
    for (int j = 0; j < _grid.sizeY(); ++j) {
        for (int i = 0; i < _grid.sizeX(); ++i) {
            if (i == 0 || j == 0 || i == _grid.sizeX() - 1 || j == _grid.sizeY() - 1) {
                _grid.setCellType(i, j, SOLID);
            }
        }
    }
}
