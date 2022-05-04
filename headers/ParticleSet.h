#pragma once
#include "LevelSet.h"
#include "MacGrid.h"

#include <memory>
#include <openvdb/openvdb.h>
#include <vector>

class LevelSet;

class Particle {
public:
    Particle();
    Particle(openvdb::Vec3d pos, openvdb::Vec3d vel);
    ~Particle();

    inline openvdb::Vec3d pos() const { return _pos; };
    inline openvdb::Vec3d vel() const { return _vel; };


    inline void setPosition(openvdb::Vec3d pos) { _pos = pos; };
    inline void setVelocity(openvdb::Vec3d vel) { _vel = vel; };

    inline void advect(float dt, openvdb::math::Transform::Ptr i2w_transform) {
        auto ws_point = i2w_transform->indexToWorld(_pos);
        ws_point += _vel * dt;
        _pos = i2w_transform->worldToIndex(ws_point);
    }

private:
    // Probably need to control the spawning time of the particle with a time_spawned variable
    openvdb::Vec3d _pos;
    openvdb::Vec3d _vel;
};
class ParticleSet {
public:
    ParticleSet(openvdb::math::Transform::Ptr i2w_transform);
    ParticleSet();
    ~ParticleSet();
    using PosType = openvdb::Vec3R;
    using ScalarType = typename PosType::value_type;

    typedef std::vector<std::unique_ptr<Particle>>::iterator iterator;
    void addParticle(std::unique_ptr<Particle> &p);
    inline void pop_back() { particles.pop_back(); };
    void removeParticle(ParticleSet::iterator p);
    void advect(float dt);
    void advectAndEnsureOutsideObstacles(LevelSet &solid_level_set, openvdb::Vec3fGrid::Ptr cpt_grid, float dt);

    inline ParticleSet::iterator begin() { return particles.begin(); };
    inline ParticleSet::iterator end() { return particles.end(); };

    inline int size() const { return particles.size(); };

    inline std::unique_ptr<Particle> &operator[](std::size_t index) { return particles[index]; };

    // Get the world-space position of the nth particle.
    // Required by rasterizeSpheres().
    inline void getPos(size_t n, openvdb::Vec3R &xyz) const {
        xyz = i2w_transform->indexToWorld(particles[n]->pos());
    };

    inline void getRadius(size_t n, ScalarType &radius) const {
        radius = this->radius;
    }

    // Get the world-space position and radius of the nth particle.
    // Required by rasterizeSpheres().
    // inline void getPosRad(size_t n, openvdb::Vec3R& xyz, openvdb::Real& radius) const {
    //     xyz[0] = 0;
    //     xyz[1] = particles[n]->posY();
    //     xyz[2] = particles[n]->posX();
    //     radius = 5;
    // };
    // // Get the world-space position, radius and velocity of the nth particle.
    // // Required by rasterizeTrails().
    // inline void getPosRadVel(size_t n, openvdb::Vec3R& xyz, openvdb::Real& radius, openvdb::Vec3R& velocity) const {
    //     xyz[0] = 0;
    //     xyz[1] = particles[n]->posY();
    //     xyz[2] = particles[n]->posX();
    //     radius = 1;
    //     velocity[0] = 0;
    //     velocity[1] = particles[n]->velY();
    //     velocity[2] = particles[n]->velX();
    // };
    // // Get the value of the nth particle's user-defined attribute (of type @c AttributeType).
    // Required only if attribute transfer is enabled in ParticlesToLevelSet.
    // inline void getAtt(size_t n, AttributeType& att) const;

    inline void clear() { particles.clear(); };
    inline void setRadius(const float radius) { this->radius = radius; }
    inline float getRadius() { return radius; }

    openvdb::math::Transform::Ptr i2w_transform;

private:
    std::vector<std::unique_ptr<Particle>> particles;
    float radius;
};