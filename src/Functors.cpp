#include "Functors.h"
#include "FluidSimulator.h"

const float epsilon = 10e-37;

void GatherTransfer::operator()(const tbb::blocked_range<int> &i_range) const {
    auto back_accessor = vel->getAccessor();
    auto vel_weight_accessor = vel_weight->getAccessor();
    openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Iterator p_it(*p_atlas);

    for (int i = i_range.begin(); i < i_range.end(); ++i) {
        for (int j = bbox.min()[1]; j < bbox.max()[1]; ++j) {
            for (int k = bbox.min()[2]; k < bbox.max()[2]; ++k) {

                // Since the kernel function is the trilinear hat in [-1,1], we only need to check 
                // cells within one cell radius 
                p_it.worldSpaceSearchAndUpdate(
                    domain._i2w_transform->indexToWorld(
                        openvdb::BBoxd(openvdb::Vec3d(i - 1, j - 1, k - 1), openvdb::Vec3d(i + 1, j + 1, k + 1))),
                    domain.particleSet());

                // Iterate over all particles contributing to the current cell, as provided by the atlas
                for (; p_it; ++p_it) {
                    auto &p = domain.particleSet()[*p_it];

                    auto coord = openvdb::Coord(i, j, k);

                    // Calculate the weight function for the current particle
                    openvdb::Vec3d weight_vec(
                        FluidSimulator::calculate_kernel_function_staggered(p.pos() - coord.asVec3d()));

                    // Update velocity and weight grids
                    if (!weight_vec.isZero()) {
                        vel_weight_accessor.setValue(coord, vel_weight_accessor.getValue(coord) + weight_vec);
                        back_accessor.setValue(coord, back_accessor.getValue(coord) + weight_vec * p.vel());
                    }
                }
            }
        }
    }
}

GatherTransfer::GatherTransfer(FluidDomain &domain,
                               openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &p_atlas,
                               openvdb::CoordBBox &bbox, const openvdb::Vec3dGrid::Ptr &vel,
                               const openvdb::Vec3dGrid::Ptr &vel_weight):
 domain(domain),
 p_atlas(p_atlas), bbox(bbox), vel_weight(vel_weight), vel(vel) {}


void ShootingTransfer::operator()(const tbb::blocked_range<int> &p_range) {
    auto back_accessor = vel->getAccessor();
    auto vel_weight_accessor = vel_weight->getAccessor();
    for (int p = p_range.begin(); p < p_range.end(); ++p) {
        auto &particle = p_set[p];

        // Since the kernel function is the trilinear hat in [-1,1], we only need to check
        // cells within one cell radius
        auto start_coord = openvdb::Coord::round(particle.pos()).offsetBy(-1, -1, -1);
        auto end_coord = start_coord.offsetBy(3, 3, 3);
        openvdb::CoordBBox bbox(start_coord, end_coord);

        for (int i = bbox.min()[0]; i < bbox.max()[0]; ++i) {
            for (int j = bbox.min()[1]; j < bbox.max()[1]; ++j) {
                for (int k = bbox.min()[2]; k < bbox.max()[2]; ++k) {
                    auto coord = openvdb::Coord(i, j, k);

                    // Calculate the weight function for the current particle
                    openvdb::Vec3d weight_vec(
                        FluidSimulator::calculate_kernel_function_staggered(particle.pos() - coord.asVec3d()));

                    // Update velocity and weight grids
                    if (!weight_vec.isZero()) {
                        vel_weight_accessor.setValue(coord, vel_weight_accessor.getValue(coord) + weight_vec);
                        back_accessor.setValue(coord, back_accessor.getValue(coord) + weight_vec * particle.vel());
                    }
                }
            }
        }
    }
}

// Constructor for tbb parallel reduce to spawn a process from an existing one
ShootingTransfer::ShootingTransfer(ShootingTransfer &x, tbb::split):
 domain(x.domain), p_set(x.p_set), vel_weight(openvdb::Vec3DGrid::create(openvdb::Vec3d(epsilon))),
 vel(openvdb::Vec3DGrid::create(openvdb::Vec3d(0))) {}

ShootingTransfer::ShootingTransfer(FluidDomain &domain, ParticleSet &p_set):
 domain(domain), p_set(p_set), vel_weight(openvdb::Vec3DGrid::create(openvdb::Vec3d(epsilon))),
 vel(openvdb::Vec3DGrid::create(openvdb::Vec3d(0))) {
    vel_weight->setTransform(domain._i2w_transform);
    vel->setGridClass(openvdb::GRID_STAGGERED);
    vel->setTransform(domain._i2w_transform);
}

void ShootingTransfer::join(const ShootingTransfer &y) {
    // Add the contributions and weights for every cell between different threads
    openvdb::tools::compSum(*vel, *y.vel);
    openvdb::tools::compSum(*vel_weight, *y.vel_weight);
}

// THIS NEEDS REVISION
void ReseedingFunctor::operator()(openvdb::tree::IteratorRange<openvdb::MaskGrid::ValueOnCIter> &range) {
    MacGrid::StaggeredSampler vel_sampler(vel->getAccessor(), vel->transform());
    openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Iterator p_it(*p_atlas);
    auto active_mask_accessor = active_mask->getAccessor();
    for (; range; ++range) {
        auto &it = range.iterator();
        auto cell_bbox = openvdb::BBoxd(it.getCoord().asVec3d() - openvdb::Vec3d(0.5, 0.5, 0.5),
                                        it.getCoord().asVec3d() + openvdb::Vec3d(0.5, 0.5, 0.5));

        p_it.worldSpaceSearchAndUpdate(domain._i2w_transform->indexToWorld(cell_bbox), domain.particleSet());

        openvdb::BBoxd particles_bbox = openvdb::BBoxd(openvdb::Vec3d(1, 1, 1), openvdb::Vec3d(-1, -1, -1));
        std::vector<int> particles_inside;
        int p_counter = 0;
        int test_counter = 0;
        for (; p_it; ++p_it) {
            auto &p = domain.particleSet()[*p_it];
            test_counter++;

            if (cell_bbox.isInside(p.pos())) {
                particles_inside.push_back(*p_it);
                p_counter++;
                // Initialise bbox
                if (particles_bbox.max() < particles_bbox.min()) {
                    particles_bbox = openvdb::BBoxd(p.pos(), p.pos());
                } else {
                    particles_bbox.expand(p.pos());
                }
            }
        }
        if (p_counter > max_particles) {
            for (auto it = particles_inside.begin(); p_counter > max_particles && it != particles_inside.end(); ++it) {
                auto particles_bbox_min = particles_bbox.min();
                auto particles_bbox_max = particles_bbox.max();
                auto particle_pos = domain.particleSet()[*it].pos();
                if (particles_bbox_min[0] < particle_pos[0] && particles_bbox_min[1] < particle_pos[1]
                    && particles_bbox_min[2] < particle_pos[2] && particle_pos[0] < particles_bbox_max[0]
                    && particle_pos[1] < particles_bbox_max[1] && particle_pos[2] < particles_bbox_max[2]) {
                    particles_to_be_deleted.push_back(*it);
                    p_counter--;
                }
            }
        } else if (p_counter < min_particles) {
            if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(-1, 0, 0))) {
                auto new_extent = cell_bbox.min();
                new_extent[0] = particles_bbox.min()[0];
                cell_bbox = openvdb::BBoxd(new_extent, cell_bbox.max());
            }
            if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(0, -1, 0))) {
                auto new_extent = cell_bbox.min();
                new_extent[1] = particles_bbox.min()[1];
                cell_bbox = openvdb::BBoxd(new_extent, cell_bbox.max());
            }
            if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(0, 0, -1))) {
                auto new_extent = cell_bbox.min();
                new_extent[2] = particles_bbox.min()[2];
                cell_bbox = openvdb::BBoxd(new_extent, cell_bbox.max());
            }
            if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(1, 0, 0))) {
                auto new_extent = cell_bbox.max();
                new_extent[0] = particles_bbox.max()[0];
                cell_bbox = openvdb::BBoxd(cell_bbox.min(), new_extent);
            }
            if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(0, 1, 0))) {
                auto new_extent = cell_bbox.max();
                new_extent[1] = particles_bbox.max()[1];
                cell_bbox = openvdb::BBoxd(cell_bbox.min(), new_extent);
            }
            if (!active_mask_accessor.isValueOn(it.getCoord().offsetBy(0, 0, 1))) {
                auto new_extent = cell_bbox.max();
                new_extent[2] = particles_bbox.max()[2];
                cell_bbox = openvdb::BBoxd(cell_bbox.min(), new_extent);
            }

            while (p_counter++ < min_particles) {
                openvdb::Vec3d new_pos =
                    cell_bbox.min()
                    + openvdb::Vec3d(dis(gen), dis(gen), dis(gen)) * (cell_bbox.max() - cell_bbox.min());
                openvdb::Vec3d new_vel = domain.grid().velInterpolatedI(vel_sampler, new_pos);
                Particle p(new_pos, new_vel, 1.01f * domain.voxelSize() * sqrt(3) / 2.f, 0);
                particles_to_be_added.push_back(p);
            }
        }
    }
}

ReseedingFunctor::ReseedingFunctor(ReseedingFunctor &x, tbb::split):
 domain(x.domain), p_atlas(x.p_atlas), active_mask(x.active_mask), vel(x.vel), min_particles(x.min_particles),
 max_particles(x.max_particles) {
    gen = std::mt19937(rd());
    dis = std::uniform_real_distribution(0., 1.);
}

ReseedingFunctor::ReseedingFunctor(FluidDomain &domain,
                                   openvdb::tools::ParticleAtlas<openvdb::tools::PointIndexGrid>::Ptr &p_atlas,
                                   openvdb::MaskGrid::Ptr active_mask, openvdb::Vec3dGrid::Ptr vel,
                                   const int min_particles, const int max_particles):
 domain(domain),
 p_atlas(p_atlas), active_mask(active_mask), vel(vel), min_particles(min_particles), max_particles(max_particles) {
    gen = std::mt19937(rd());
    dis = std::uniform_real_distribution(0., 1.);
}

void ReseedingFunctor::join(const ReseedingFunctor &y) {
    particles_to_be_added.reserve(particles_to_be_added.size() + y.particles_to_be_added.size());
    particles_to_be_deleted.reserve(particles_to_be_deleted.size() + y.particles_to_be_deleted.size());

    particles_to_be_added.insert(particles_to_be_added.end(), y.particles_to_be_added.begin(),
                                 y.particles_to_be_added.end());
    particles_to_be_deleted.insert(particles_to_be_deleted.end(), y.particles_to_be_deleted.begin(),
                                   y.particles_to_be_deleted.end());
}