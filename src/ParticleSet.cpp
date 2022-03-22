#include "ParticleSet.h"

Particle::Particle() {}
Particle::Particle(float pos_x, float pos_y, float vel_x, float vel_y, float temperature, float concentration):
 pos_x(pos_x), pos_y(pos_y), vel_x(vel_x), vel_y(vel_y), _temperature(temperature), _concentration(concentration) {};
Particle::~Particle() {}

ParticleSet::ParticleSet() {};
ParticleSet::~ParticleSet() {};

void ParticleSet::addParticle(std::unique_ptr<Particle> &p) { particles.push_back(std::move(p)); }
void ParticleSet::removeParticle(ParticleSet::iterator p) {
    particles.erase(p);
}
void ParticleSet::advect(float dt) {
    for (auto it = particles.begin(); it != particles.end(); ++it) (*it)->advect(dt);
}
// Needs to change
void ParticleSet::advectAndEnsureOutsideObstacles(MacGrid &grid, float dt) {
    for (auto it = particles.begin(); it != particles.end(); it++) {
        (*it)->advect(dt);
        int i = (*it)->posX() / grid.deltaX();
        int j = (*it)->posY() / grid.deltaY();
        if (grid.cellType(i, j) == SOLID) (*it)->advect(-dt);
    }
}