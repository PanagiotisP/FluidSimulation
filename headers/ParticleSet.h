#pragma once
#include "LevelSet.h"
#include "MacGrid.h"

#include <memory>
#include <openvdb/openvdb.h>
#include <openvdb/tools/VelocityFields.h>
#include <vector>

class LevelSet;

class Particle {
public:
    Particle();
    Particle(openvdb::Vec3d pos, openvdb::Vec3d vel, float radius, float mass);
    ~Particle();

    // Getters
    inline openvdb::Vec3d pos() const { return _pos; };
    inline openvdb::Vec3d vel() const { return _vel; };
    inline float radius() const { return _radius; };
    inline float mass() const { return _mass; };

    // Setters
    inline void setPosition(openvdb::Vec3d pos) { _pos = pos; };
    inline void setVelocity(openvdb::Vec3d vel) { _vel = vel; };

    // Advect a particle's world position using an RK method, provided as function argument
    inline void advect(float dt, openvdb::math::Transform::Ptr i2w_transform,
                       const openvdb::tools::VelocityIntegrator<openvdb::Vec3dGrid, 1> &rk_integrator) {
        auto ws_point = i2w_transform->indexToWorld(_pos);
        rk_integrator.rungeKutta<1, openvdb::Vec3d>(dt, ws_point);
        _pos = i2w_transform->worldToIndex(ws_point);
    }

private:
    openvdb::Vec3d _pos;
    openvdb::Vec3d _vel;
    float _radius;
    float _mass;
};

class ParticleSet {
public:
    ParticleSet();
    ParticleSet(openvdb::math::Transform::Ptr i2w_transform);
    ~ParticleSet();

    using PosType = openvdb::Vec3R;
    using ScalarType = typename PosType::value_type;

    typedef std::vector<Particle>::iterator iterator;

    // Advect particles and project those which end up inside solids on the closest point of the solid's surface
    void advectAndEnsureOutsideObstacles(openvdb::Vec3dGrid::Ptr vel_grid, openvdb::Vec3fGrid::Ptr cpt_grid,
                                         LevelSet &solid_level_set, float dt);

    inline void addParticle(const Particle &p) {
        _particles.push_back(p);
        _minRadius = min(_minRadius, p.radius());
    };
    inline void pop_back() { _particles.pop_back(); };
    inline void removeParticle(ParticleSet::iterator p) { _particles.erase(p); };

    inline ParticleSet::iterator begin() { return _particles.begin(); };
    inline ParticleSet::iterator end() { return _particles.end(); };

    inline int size() const { return _particles.size(); };

    inline float getMinRadius() { return _minRadius; };

    inline Particle &operator[](std::size_t index) { return _particles[index]; };

    // Get the world-space position of the nth particle.
    // Required by rasterizeSpheres().
    inline void getPos(size_t n, openvdb::Vec3R &xyz) const {
        xyz = _i2w_transform->indexToWorld(_particles[n].pos());
    };

    // Get the world-space radius of the nth particle.
    // Required by rasterizeSpheres().
    inline void getRadius(size_t n, ScalarType &radius) const {
        if (_zeroRadius) {
            radius = 0;
        } else {
            radius = _particles[n].radius();
        }
    }

    // Get the world-space position and radius of the nth particle.
    // Required by rasterizeSpheres().
    inline void getPosRad(size_t n, openvdb::Vec3R &xyz, openvdb::Real &radius) const {
        xyz = _i2w_transform->indexToWorld(_particles[n].pos());
        if (_zeroRadius) {
            radius = 0;
        } else {
            radius = _particles[n].radius();
        }
    };

    // Empty set of all particles
    inline void clear() { _particles.clear(); };

    inline void setZeroRadius(bool zeroRadius) { _zeroRadius = zeroRadius; };
    openvdb::math::Transform::Ptr _i2w_transform;

private:
    std::vector<Particle> _particles;

    // This is only used for ParticlesToLevelSet RMin parameter,
    // so as to not have particles that get ignored due to small radius.
    // We don't really need to increase it (e.g. on particleRemove).
    // We just want, at any point in time, to be less than or equal to the min particles' radius
    float _minRadius = INT_MAX;
    bool _zeroRadius = false;
};