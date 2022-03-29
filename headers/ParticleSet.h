#pragma once
#include "LevelSet.h"
#include "MacGrid.h"

#include <memory>
#include <vector>

class Particle {
public:
    Particle();
    Particle(float pos_x, float pos_y, float vel_x, float vel_y, float temperature, float concentration);
    ~Particle();

    inline float posX() const { return pos_x; };
    inline float posY() const { return pos_y; };
    inline float velX() const { return vel_x; };
    inline float velY() const { return vel_y; };
    inline float temperature() const { return _temperature; };
    inline float concentration() const { return _concentration; };


    inline void setPosition(float x, float y) {
        pos_x = x;
        pos_y = y;
    };
    inline void setVelocity(float vel_x, float vel_y) {
        this->vel_x = vel_x;
        this->vel_y = vel_y;
    };
    inline void setTemperature(float val) { _temperature = val; };
    inline void setConcentration(float val) { _concentration = val; };

    inline void advect(float dt) {
        pos_x += vel_x * dt;
        pos_y += vel_y * dt;
    }

private:
// Probably need to control the spawning time of the particle with a time_spawned variable
    float pos_x, pos_y;
    float vel_x, vel_y;
    float _temperature;
    float _concentration;
};
class ParticleSet {
public:
    ParticleSet();
    ~ParticleSet();

    typedef std::vector<std::unique_ptr<Particle>>::iterator iterator;
    void addParticle(std::unique_ptr<Particle> &p);
    inline void pop_back() { particles.pop_back(); };
    void removeParticle(ParticleSet::iterator p);
    void advect(float dt);
    void advectAndEnsureOutsideObstacles(LevelSet& solid_level_set, MacGrid &mac_grid, float dt);

    inline ParticleSet::iterator begin() { return particles.begin(); };
    inline ParticleSet::iterator end() { return particles.end(); };

    inline int size() const { return particles.size(); };
    inline void clear() { particles.clear(); };

private:
    std::vector<std::unique_ptr<Particle>> particles;
};