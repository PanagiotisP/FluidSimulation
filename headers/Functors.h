#include "FluidDomain.h"
#include "ParticleSet.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/ParticleAtlas.h>
#include <tbb/parallel_reduce.h>

// Class for Gather Particles to Grid Transfer
// In this case, we iterate over every cell and gather (so the name) the 
// contribution of all particles inside it.
class GatherTransfer {
public:
    // These two grids have preallocated topology of the final grid,
    // meaning that we can perform parallel writes with the ValueAccessor
    const openvdb::Vec3dGrid::Ptr &vel_weight;
    const openvdb::Vec3dGrid::Ptr &vel;
    
    // ParticleAtlas is an accelerator structure that finds the particles within a radius efficiently
    openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &p_atlas;
    openvdb::CoordBBox &bbox;
    FluidDomain &domain;

    void operator()(const tbb::blocked_range<int> &i_range) const;

    GatherTransfer(FluidDomain &domain, openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &p_atlas,
                   openvdb::CoordBBox &bbox, const openvdb::Vec3dGrid::Ptr &vel,
                   const openvdb::Vec3dGrid::Ptr &vel_weight);
};

// Class for Shooting Particles to Grid Transfer
// In this case, we iterate over every particle and distribute (shoot) its 
// contribution to all cells that are affected by it, as determined by the kernel function.
class ShootingTransfer {
public:
    // These two grids store the results of some particles and are reduced-combined together at the end
    openvdb::Vec3dGrid::Ptr vel_weight;
    openvdb::Vec3dGrid::Ptr vel;
    ParticleSet &p_set;
    FluidDomain &domain;

    void operator()(const tbb::blocked_range<int> &p_range);
    
    // Constructor for tbb parallel reduce to spawn a process from an existing one
    ShootingTransfer(ShootingTransfer &x, tbb::split);

    ShootingTransfer(FluidDomain &domain, ParticleSet &p_set);

    // Join method for parallel reduce, responsible to perform the reduce operation
    void join(const ShootingTransfer &y);
};

class ReseedingFunctor {
    std::random_device rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen;      // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis;
    const int max_particles;
    const int min_particles;
    openvdb::MaskGrid::Ptr active_mask;
    openvdb::Vec3dGrid::Ptr vel;
    FluidDomain &domain;

public:
    std::vector<int> particles_to_be_deleted;
    std::vector<Particle> particles_to_be_added;

    openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &p_atlas;
    
    void operator()(openvdb::tree::IteratorRange<openvdb::MaskGrid::ValueOnCIter> &range);

    ReseedingFunctor(ReseedingFunctor &x, tbb::split);

    ReseedingFunctor(FluidDomain &domain, openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &p_atlas,
                     openvdb::MaskGrid::Ptr active_mask, openvdb::Vec3dGrid::Ptr vel, const int min_particles,
                     const int max_particles);

    void join(const ReseedingFunctor &y);
};

// Class for Density Calculation
// In this case, we iterate over every cell and gather (so the name) the 
// contribution of all particles inside it, while estimating the volume from face fractions
class DensityCalculator {
public:
    // These grid have preallocated topology of the final grid,
    // meaning that we can perform parallel writes with the ValueAccessor
    const openvdb::FloatGrid::Ptr &_density_grid;

    // ParticleAtlas is an accelerator structure that finds the particles within a radius efficiently
    openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &_p_atlas;
    openvdb::CoordBBox &_bbox;
    FluidDomain &_domain;

    void operator()(openvdb::tree::IteratorRange<openvdb::MaskGrid::ValueOnCIter>) const;

    DensityCalculator(FluidDomain &domain, openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &p_atlas,
                      openvdb::CoordBBox &bbox, const openvdb::FloatGrid::Ptr &density_grid);
};
