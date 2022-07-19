#include "FluidSource.h"

FluidSource::FluidSource():
 _spawning_region(openvdb::createLevelSet<openvdb::FloatGrid>()), _max_spawning_duration(0) {}
FluidSource::FluidSource(LevelSet spawning_region, openvdb::Vec3d vel, int particle_generation_rate,
                         LevelSet &solidLevelSet, float max_spawning_duration):
 _vel(vel),
 _spawning_region(spawning_region), _particle_generation_rate(particle_generation_rate),
 _max_spawning_duration(max_spawning_duration), _spawning_duration_left(max_spawning_duration) {
    this->_spawning_region.getLevelSet()->tree().combineExtended(solidLevelSet.getLevelSet()->deepCopy()->tree(),
                                                                 IntersectSolidFluidLevelSets::intersect, true);
}

FluidSource::~FluidSource() {}

void FluidSource::update(MacGrid &grid, ParticleSet &particle_set, float dt) {
    if (_spawning_duration_left > 0) {
        // Spawn particles based on the rate provided
        // To handle non integer amount of to-be-spawned particles
        // a particle is spawned with a probability of the decimal part of the amount

        // If spawning_duration_left is less that dt, use that instead to find spawned particles amount
        float spawn_amount = min(dt, _spawning_duration_left) * _particle_generation_rate;
        int spawn_counter = static_cast<int>(spawn_amount);
        float probability_threshold = spawn_amount - static_cast<float>(spawn_counter);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0., 1.);

        auto coordBBox = _spawning_region.getActiveCoordBBox();
        auto accessor = _spawning_region.getAccessor();
        std::shared_ptr<LevelSet::Sampler> sampler(_spawning_region.getSampler(_spawning_region.getAccessor()));
        // LevelSet::Sampler sampler(accessor, _spawning_region.getLevelSet()->transform());

        MacGrid::Sampler vel_sampler(grid.velFront()->getAccessor(), _spawning_region.getLevelSet()->transform());

        assert(coordBBox.min() < coordBBox.max());

        // Particles are spawned for each cell (partially) covered by the spawning region,
        // as determined by the sdf value.
        // Initial velocity is interpolated from the grid and the position is set as the cell's center,
        // offseted by a random vector to, in order to uniformly sample the whole cell.
        for (auto it = coordBBox.beginZYX(); it != coordBBox.endZYX(); ++it) {
            auto cellCenterCoord = (*it);
            if (accessor.getValue(cellCenterCoord) < 0) {
                int temp_spawn_counter = spawn_counter;
                if (dis(gen) < probability_threshold) temp_spawn_counter++;

                // Include a max tries for the spawning of a particle to avoid endless loops
                int tries = 0;
                while (temp_spawn_counter > 0 && tries < 8 * temp_spawn_counter) {
                    tries++;
                    openvdb::Vec3d new_pos = cellCenterCoord.asVec3d() + openvdb::Vec3d(dis(gen), dis(gen), dis(gen));

                    // Spawn only inside fluid cells
                    if (_spawning_region.valueInterpolatedI(sampler, new_pos) < 0) {
                        // TODO Need to check this addition of spawning and interpolated velocity
                        openvdb::Vec3d new_vel = grid.velInterpolatedI(vel_sampler, new_pos) + _vel;
                        // As for now, radius is the same for all particles. To be changed later.
                        Particle p(new_pos, new_vel, 1.01f * grid.voxelSize() * sqrt(3) / 2.f, 0);
                        particle_set.addParticle(p);
                        temp_spawn_counter--;
                    }
                }
            }
        }
        _spawning_duration_left -= dt;
    }
}