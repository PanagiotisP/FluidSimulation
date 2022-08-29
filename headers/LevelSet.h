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

using DistanceFunction = std::function<float(float, float)>;

class ParticleSet;

class LevelSet {
public:
    using Sampler = openvdb::tools::GridSampler< openvdb::FloatGrid::Accessor, openvdb::tools::BoxSampler >;
    using Accessor = openvdb::FloatGrid::Accessor;

    ~LevelSet();
    LevelSet(openvdb::FloatGrid::Ptr level_set);

    // Make a level set by taking the union-of-balls from the particles' data
    // This means that each particle is treated as a ball of a given radius
    // and the level set is the union of these balls
    void unionOfBalls(ParticleSet &points, float voxel_size);

    // Value interpolation in index space
    float valueInterpolatedI(Sampler &sampler, float x, float y, float z);
    float valueInterpolatedI(std::shared_ptr<Sampler> &sampler, openvdb::Vec3R point);

    // Value interpolation in world space
    float valueInterpolatedW(Sampler &sampler, float x, float y, float z);
    float valueInterpolatedW(Sampler &sampler, openvdb::Vec3R point);


    // CSG Union. The union of self and argument level_sets is stored in self level_set
    void unionLevelSetInPlace(LevelSet &level_set);
    // CSG Intersection. The intersection of self and argument level_sets is stored in self level_set
    void intersectionLevelSetInPlace(LevelSet &level_set);
    // CSG Difference. The difference of self and argument level_sets is stored in self level_set
    void differenceLevelSetInPlace(LevelSet &level_set);

    // CSG Union. Returns a new LevelSet with the result of the operation.
    friend LevelSet unionLevelSet(const LevelSet &level_set_a, const LevelSet &level_set_b);
    // CSG Intersection. Returns a new LevelSet with the result of the operation.
    friend LevelSet intersectionLevelSet(const LevelSet &level_set_a, const LevelSet &level_set_b);
    // CSG Difference. Returns a new LevelSet with the result of the operation.
    friend LevelSet differenceLevelSet(const LevelSet &level_set_a, const LevelSet &level_set_b);

    // Advect level set
    void advect(openvdb::Vec3dGrid::Ptr vel, float dt);

    // Change the sign of each level set entry
    // For sdf levelsets this is equivalent to turning the relative shape inside out
    void invert();

    inline openvdb::CoordBBox getActiveCoordBBox() { return _level_set->evalActiveVoxelBoundingBox(); };

    inline Accessor &getAccessor() {
        _accessor = _level_set->getAccessor();
        return _accessor;
    };

    inline const std::shared_ptr<Sampler> getSampler(const Accessor &accessor) {
        return std::make_shared<Sampler>(accessor, _level_set->transform());
    };

    inline openvdb::FloatGrid::Ptr getLevelSet() const { return _level_set; };
    inline void setLevelSet(openvdb::FloatGrid::Ptr const &level_set) { _level_set = level_set; };

    friend std::ostream &operator<<(std::ostream &out, LevelSet &l);

private:
    openvdb::FloatGrid::Ptr _level_set;
    Accessor _accessor;
    Sampler _sampler;
};

struct IntersectSolidFluidLevelSets {
    static inline void intersect(openvdb::CombineArgs<float> &args) {
        // Args a fluid, Args b solid
        // If inside the solid then either the fluid is not touching (max = fluid_val), or it is touching and is
        // -solid distance away
        float result;
        if (args.b() <= 0)
            result = max(args.a(), -args.b());
        else
            result = args.a();
        args.setResult(result);
        args.setResultIsActive(args.aIsActive());
    }
};

#endif