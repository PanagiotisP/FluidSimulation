#ifndef LEVELSET_H
#define LEVELSET_H

#include "ParticleSet.h"
#include "util.h"
#include "vec.h"

#include <functional>
#include <iostream>
#include <math.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

typedef std::function<float(float, float)> DistanceFunction;

class ParticleSet;

class LevelSet {
public:
    typedef openvdb::tools::GridSampler< openvdb::FloatGrid::Accessor, openvdb::tools::BoxSampler > Sampler;
    typedef openvdb::FloatGrid::Accessor FloatAccessor;
    LevelSet(openvdb::FloatGrid::Ptr level_set);
    ~LevelSet();

    void unionOfBalls(ParticleSet &points, float voxel_size);

    inline openvdb::FloatGrid::Accessor getAccessor() { return _level_set->getAccessor(); };
    inline Sampler getBoxSampler() { return Sampler(_level_set->getAccessor(), _level_set->transform()); };
    float valueInterpolatedI(Sampler &sampler, float x, float y, float z);
    float valueInterpolatedI(Sampler &sampler, openvdb::Vec3R point);
    float valueInterpolatedW(Sampler &sampler, float x, float y, float z);
    float valueInterpolatedW(Sampler &sampler, openvdb::Vec3R point);

    inline openvdb::CoordBBox getActiveCoordBBox() { return _level_set->evalActiveVoxelBoundingBox(); };

    void unionLevelSet(LevelSet &level_set);
    void intersectionLevelSet(LevelSet &level_set);
    void differenceLevelSet(LevelSet &level_set);

    inline openvdb::FloatGrid::Ptr getLevelSet() { return _level_set; };
    inline void setLevelSet(openvdb::FloatGrid::Ptr const & level_set) { _level_set = level_set; };

    // Takes an inside-outside level set that is far from sdf and produces sdf
    void redistance();

    LevelSet invert();
    friend std::ostream &operator<<(std::ostream &out, LevelSet &l);
    static double fraction_inside(double phi_left, double phi_right);
    static double fraction_inside(double phi_bl, double phi_br, double phi_tl, double phi_tr);
    static double fraction_inside(double phi_bl, double phi_br, double phi_tl, double phi_tr, double phi_center);
    static void cycle_array(double *arr, int size);

private:
    openvdb::FloatGrid::Ptr _level_set;
};

#endif