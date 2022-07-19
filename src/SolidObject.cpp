#pragma once
#include "SolidObject.h"

#include "LevelSet.h"
#include "openvdb/tools/FastSweeping.h"
#include "openvdb/tools/GridTransformer.h"
#include "openvdb/tools/LevelSetSphere.h"
#include "openvdb/tools/VelocityFields.h"


//  SPHERE FOR NOW
SolidObject::SolidObject(float radius, float voxel_size, double density, openvdb::math::Transform::Ptr xform,
                         openvdb::Vec3d center, openvdb::Vec3d vel_translational, openvdb::Vec3d vel_angular):
 _center_of_mass(center / voxel_size),
 xform(xform), _vel_translational(vel_translational), _vel_angular(vel_angular), _density(density),
 _level_set(openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(radius, center, voxel_size, 5)),
 _u_weight_accessor(_u_weights->getAccessor()), _v_weight_accessor(_v_weights->getAccessor()),
 _w_weight_accessor(_w_weights->getAccessor()) {
    _level_set.getLevelSet()->setTransform(xform);
    _u_weights->setTransform(xform);
    _v_weights->setTransform(xform);
    _w_weights->setTransform(xform);

    float mass = 4 / 3 * 3.14 * pow(radius, 3) * _density;

    Eigen::Matrix<double, 6, 6> mass_matrix;

    mass_matrix.setZero();

    mass_matrix(0, 0) = mass;
    mass_matrix(1, 1) = mass;
    mass_matrix(2, 2) = mass;

    // Inertia tensor
    mass_matrix(3, 3) = 2 * mass * sqr(radius) / 5;
    mass_matrix(4, 4) = 2 * mass * sqr(radius) / 5;
    mass_matrix(5, 5) = 2 * mass * sqr(radius) / 5;

    _mass_matrix_inverse = mass_matrix.inverse().sparseView();
}

Eigen::Vector<double, 6> SolidObject::combinedVelocityVector() {
    Eigen::Vector<double, 6> vec;
    vec[0] = _vel_translational[0];
    vec[1] = _vel_translational[1];
    vec[2] = _vel_translational[2];

    vec[3] = _vel_angular[0];
    vec[4] = _vel_angular[1];
    vec[5] = _vel_angular[2];

    return vec;
}

void SolidObject::setVelocityVector(Eigen::Vector<double, 6> vel_vector) {
    _vel_translational[0] = vel_vector[0];
    _vel_translational[1] = vel_vector[1];
    _vel_translational[2] = vel_vector[2];

    _vel_angular[0] = vel_vector[3];
    _vel_angular[1] = vel_vector[4];
    _vel_angular[2] = vel_vector[5];
}

void SolidObject::generateVelocityField() {
    vel_field = openvdb::createGrid<openvdb::Vec3dGrid>(_vel_translational);
    vel_field->setTransform(xform);
    auto voxel_size = _level_set.getLevelSet()->voxelSize()[0];
    auto coordBBox = _level_set.getActiveCoordBBox();
    auto vel_accessor = vel_field->getAccessor();
    if (coordBBox.min() < coordBBox.max()) {
        for (auto it = coordBBox.beginXYZ(); it != coordBBox.endXYZ(); ++it) {
            openvdb::Vec3d staggered_velocity;
            staggered_velocity[0] = _vel_angular.cross(voxel_size * ((*it).asVec3d() - openvdb::Vec3d(0.5, 0, 0)))[0];
            staggered_velocity[1] = _vel_angular.cross(voxel_size * ((*it).asVec3d() - openvdb::Vec3d(0, 0.5, 0)))[1];
            staggered_velocity[2] = _vel_angular.cross(voxel_size * ((*it).asVec3d() - openvdb::Vec3d(0, 0, 0.5)))[2];

            vel_accessor.setValue((*it), _vel_translational + staggered_velocity);
        }
    }
}

void SolidObject::update_positions_orientations(float dt) {
    auto voxel_size = _level_set.getLevelSet()->voxelSize()[0];
    // Point is affected only by translation velocity
    openvdb::Vec3d scaled_com = _center_of_mass * voxel_size;
    scaled_com += dt * _vel_translational;
    _center_of_mass = scaled_com / voxel_size;

    // for (int i = coordBBox.min()[0]; i < coordBBox.max()[0]; ++i) {
    //     std::cout << "\n" << i << std::endl;
    //     for (int j = coordBBox.max()[1]; j >= coordBBox.min()[1]; --j) {
    //         std::cout << std::endl;
    //         for (int k = coordBBox.min()[2]; k < coordBBox.max()[2]; ++k) {
    //             std::cout << vel_field->getAccessor().getValue(openvdb::Coord(i, j, k)) << " ";
    //         }
    //     }
    // }

    _level_set.advect(vel_field, dt);
}

void SolidObject::compute_face_fractions() {

    // Resample solid level set to get values at bottom left corner of voxels
    auto resampledSolidLevelSet = openvdb::createLevelSet<openvdb::FloatGrid>(_level_set.getLevelSet()->voxelSize()[0]);
    openvdb::tools::GridTransformer transformer(
        _level_set.getLevelSet()->transformPtr()->baseMap()->getAffineMap()->getMat4()
        * resampledSolidLevelSet->transformPtr()->baseMap()->getAffineMap()->getMat4().inverse());

    transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(*_level_set.getLevelSet(),
                                                                              *resampledSolidLevelSet);
    auto border_mask = openvdb::tools::extractIsosurfaceMask(*_level_set.getLevelSet(), 0);
    auto inside_mask = openvdb::tools::extractEnclosedRegion(*_level_set.getLevelSet(), 0);
    inside_mask->tree().topologyUnion(border_mask->tree());

    openvdb::tools::dilateActiveValues(inside_mask->tree(), 1, openvdb::tools::NN_FACE_EDGE_VERTEX);

    _u_weights->tree().topologyUnion(inside_mask->tree());
    _v_weights->tree().topologyUnion(inside_mask->tree());
    _w_weights->tree().topologyUnion(inside_mask->tree());

    _u_weights->tree().voxelizeActiveTiles();
    _v_weights->tree().voxelizeActiveTiles();
    _w_weights->tree().voxelizeActiveTiles();

    openvdb::tree::IteratorRange<openvdb::BoolGrid::ValueOnCIter> range(inside_mask->cbeginValueOn());
    tbb::parallel_for(range, [&](openvdb::tree::IteratorRange<openvdb::BoolGrid::ValueOnCIter> range) {
        auto targetAccessor = resampledSolidLevelSet->getAccessor();
        auto u_weights_accessor = _u_weights->getAccessor();
        auto v_weights_accessor = _v_weights->getAccessor();
        auto w_weights_accessor = _w_weights->getAccessor();
        for (; range; ++range) {
            auto coord = range.iterator().getCoord();
            auto point_000 = targetAccessor.getValue(coord.offsetBy(0, 0, 0));
            auto point_001 = targetAccessor.getValue(coord.offsetBy(0, 0, 1));
            auto point_010 = targetAccessor.getValue(coord.offsetBy(0, 1, 0));
            auto point_011 = targetAccessor.getValue(coord.offsetBy(0, 1, 1));
            auto point_100 = targetAccessor.getValue(coord.offsetBy(1, 0, 0));
            auto point_101 = targetAccessor.getValue(coord.offsetBy(1, 0, 1));
            auto point_110 = targetAccessor.getValue(coord.offsetBy(1, 1, 0));

            u_weights_accessor.setValue(coord, fraction_inside(point_000, point_010, point_001, point_011));
            v_weights_accessor.setValue(coord, fraction_inside(point_000, point_001, point_100, point_101));
            w_weights_accessor.setValue(coord, fraction_inside(point_000, point_010, point_100, point_110));
        }
    });
}

void SolidObject::add_forces(float dt) {
    // Only gravity
    _vel_translational[1] -= 9.81f * dt;
}

void SolidObject::constructJMatrix(int idx, openvdb::Coord coord, float voxelSize) {
    _j_matrix.insert(0, idx) =
        sqr(voxelSize) * (_u_weight_accessor.getValue(coord.offsetBy(1, 0, 0)) - _u_weight_accessor.getValue(coord));
    _j_matrix.insert(1, idx) =
        sqr(voxelSize) * (_v_weight_accessor.getValue(coord.offsetBy(0, 1, 0)) - _v_weight_accessor.getValue(coord));
    _j_matrix.insert(2, idx) =
        sqr(voxelSize) * (_w_weight_accessor.getValue(coord.offsetBy(0, 0, 1)) - _w_weight_accessor.getValue(coord));

    _j_matrix.insert(3, idx) =
        sqr(voxelSize)
        * (-voxelSize * (coord[2] - _center_of_mass[2])
               * (_v_weight_accessor.getValue(coord.offsetBy(0, 1, 0)) - _v_weight_accessor.getValue(coord))
           + voxelSize * (coord[1] - _center_of_mass[1])
                 * (_w_weight_accessor.getValue(coord.offsetBy(0, 0, 1)) - _w_weight_accessor.getValue(coord)));
    _j_matrix.insert(4, idx) =
        sqr(voxelSize)
        * (-voxelSize * (coord[0] - _center_of_mass[0])
               * (_w_weight_accessor.getValue(coord.offsetBy(0, 0, 1)) - _w_weight_accessor.getValue(coord))
           + voxelSize * (coord[2] - _center_of_mass[2])
                 * (_u_weight_accessor.getValue(coord.offsetBy(1, 0, 0)) - _u_weight_accessor.getValue(coord)));
    _j_matrix.insert(5, idx) =
        sqr(voxelSize)
        * (-voxelSize * (coord[1] - _center_of_mass[1])
               * (_u_weight_accessor.getValue(coord.offsetBy(1, 0, 0)) - _u_weight_accessor.getValue(coord))
           + voxelSize * (coord[0] - _center_of_mass[0])
                 * (_v_weight_accessor.getValue(coord.offsetBy(0, 1, 0)) - _v_weight_accessor.getValue(coord)));
}