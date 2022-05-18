#pragma once

#ifndef MAC_GRID_H
    #define MAC_GRID_H

//#include <glm/glm.hpp>

    #include "util.h"

    #include <assert.h>
    #include <iostream>
    #include <memory>
    #include <openvdb/openvdb.h>
    #include <openvdb/tools/Interpolation.h>
// Some comments about the MAC Grid implementation
// The int arguments in member functions refer to cell indices of the initial MacGrid.
// As it is staggered, this means that for scalar values taken at the middle of the cell the value (i,j) is returned,
// while for values taken at half indices an interpolation of two adjacent cells is taken
// The float values refer to world coordinates.
// In interpolating functions that get world coordinates as input the following rules apply
// For interpolating the x-component of velocity, the world input should have an offset of (0,0.5), so that makes it for
// the center of a grid cell (i,j) (i*dx, (j+0.5)*dy) For interpolating the y-component of velocity, the world input
// should have an offset of (0.5,0), so that makes it for the center of a grid cell (i,j) ((i+0.5)*dx, j*dy) For
// interpolating scalar values taken at the center of the cell (i,j), the world input should have an offset of
// (0.5,0.5), so that makes ((i+0.5)*dx, (j+0.5)*dy) (the last one is nothing more than getting the correct coordinates
// of the center of a cell)


class MacGrid {
public:
    MacGrid(openvdb::math::Transform::Ptr i2w_transform);
    ~MacGrid();

    typedef openvdb::Vec3dGrid::Accessor Accessor;
    typedef openvdb::tools::GridSampler<openvdb::Vec3dGrid::Accessor, openvdb::tools::StaggeredBoxSampler> Sampler;

    inline openvdb::Vec3dGrid::Ptr velFront() { return _vel_front; }
    inline openvdb::Vec3dGrid::Ptr velBack() { return _vel_back; }
    inline openvdb::Vec3dGrid::Ptr velPrev() { return _vel_prev; }
    inline openvdb::Vec3dGrid::Ptr velDiff() { return _vel_diff; }
    inline openvdb::MaskGrid::Ptr validUFront() { return valid_mask_u_front_buffer; };
    inline openvdb::MaskGrid::Ptr validUBack() { return valid_mask_u_back_buffer; };
    inline openvdb::MaskGrid::Ptr validVFront() { return valid_mask_v_front_buffer; };
    inline openvdb::MaskGrid::Ptr validVBack() { return valid_mask_v_back_buffer; };
    inline openvdb::MaskGrid::Ptr validWFront() { return valid_mask_w_front_buffer; };
    inline openvdb::MaskGrid::Ptr validWBack() { return valid_mask_w_back_buffer; };

    inline openvdb::FloatGrid::Ptr uWeights() { return u_weights; };
    inline openvdb::FloatGrid::Ptr vWeights() { return v_weights; };
    inline openvdb::FloatGrid::Ptr wWeights() { return w_weights; };
    

    inline void setValidUBack(openvdb::MaskGrid::Ptr mask) { valid_mask_u_back_buffer = mask; }
    inline void setValidWBack(openvdb::MaskGrid::Ptr mask) { valid_mask_w_back_buffer = mask; }
    inline void setValidVBack(openvdb::MaskGrid::Ptr mask) { valid_mask_v_back_buffer = mask; }
    inline void setVelBack(openvdb::Vec3dGrid::Ptr vel_back) { _vel_back = vel_back; }

    void updatePreviousVelocityBuffer();
    void updateDiffBuffers();

    // Getters
    inline openvdb::Vec3d velHalfIndexed(Accessor &accessor, openvdb::Coord coord) const { return accessor.getValue(coord); };
    inline openvdb::Vec3d velHalfIndexed(Accessor &accessor, int k, int j, int i) const {
        return accessor.getValue(openvdb::Coord(k, j, i));
    };

    inline openvdb::Vec3d velInterpolatedI(Sampler &sampler, openvdb::Vec3d isPoint) const {
        return sampler.isSample(isPoint);
    };
    inline openvdb::Vec3d velInterpolatedI(Sampler &sampler, double x, double y, double z) const {
        return sampler.isSample(openvdb::Vec3d(x, y, z));
    };

    inline openvdb::Vec3d velInterpolatedW(Sampler &sampler, openvdb::Vec3d wsPoint) const {
        return sampler.wsSample(wsPoint);
    };

    inline openvdb::Vec3d velInterpolatedW(Sampler &sampler, double x, double y, double z) const {
        return sampler.wsSample(openvdb::Vec3d(x, y, z));
    };

    // glm::dmat2 computeVelocityGradientMatrix(int i, int j);

    // Setters
    inline void setVelHalfIndexed(Accessor &accessor, openvdb::Coord coord, openvdb::Vec3d vel) { accessor.setValue(coord, vel); };

    inline void setVelHalfIndexed(Accessor &accessor, int k, int j, int i, openvdb::Vec3d vel) {
        accessor.setValue(openvdb::Coord(k, j, i), vel);
    };

    inline void setVelXHalfIndexed(Accessor &accessor, openvdb::Coord coord, double vel) {
        openvdb::Vec3d value = accessor.getValue(coord);
        value[0] = vel;
        accessor.setValue(coord, value);
    };

    inline void setVelXHalfIndexed(Accessor &accessor, int k, int j, int i, double vel) {
        openvdb::Vec3d value = accessor.getValue(openvdb::Coord(k, j, i));
        value[0] = vel;
        accessor.setValue(openvdb::Coord(k, j, i), value);
    };

    inline void setVelYHalfIndexed(Accessor &accessor, openvdb::Coord coord, double vel) {
        openvdb::Vec3d value = accessor.getValue(coord);
        value[1] = vel;
        accessor.setValue(coord, value);
    };

    inline void setVelYHalfIndexed(Accessor &accessor, int k, int j, int i, double vel) {
        openvdb::Vec3d value = accessor.getValue(openvdb::Coord(k, j, i));
        value[1] = vel;
        accessor.setValue(openvdb::Coord(k, j, i), value);
    };

    inline void setVelZHalfIndexed(Accessor &accessor, openvdb::Coord coord, double vel) {
        openvdb::Vec3d value = accessor.getValue(coord);
        value[2] = vel;
        accessor.setValue(coord, value);
    };

    inline void setVelZHalfIndexed(Accessor &accessor, int k, int j, int i, double vel) {
        openvdb::Vec3d value = accessor.getValue(openvdb::Coord(k, j, i));
        value[2] = vel;
        accessor.setValue(openvdb::Coord(k, j, i), value);
    };

    inline void setVel(Accessor &accessor, int k, int j, int i, openvdb::Vec3d vel) {
        accessor.setValue(openvdb::Coord(k, j, i), vel);
        accessor.setValue(openvdb::Coord(k, j, i + 1), vel);
        accessor.setValue(openvdb::Coord(k, j + 1, i), vel);
        accessor.setValue(openvdb::Coord(k, j + 1, i + 1), vel);
        accessor.setValue(openvdb::Coord(k + 1, j, i), vel);
        accessor.setValue(openvdb::Coord(k + 1, j, i + 1), vel);
        accessor.setValue(openvdb::Coord(k + 1, j + 1, i), vel);
        accessor.setValue(openvdb::Coord(k + 1, j + 1, i + 1), vel);
    };

    inline void setVel(Accessor &accessor, openvdb::Coord coord, openvdb::Vec3d vel) {
        accessor.setValue(coord, vel);
        accessor.setValue(coord.offsetBy(0, 0, 1), vel);
        accessor.setValue(coord.offsetBy(0, 1, 0), vel);
        accessor.setValue(coord.offsetBy(0, 1, 1), vel);
        accessor.setValue(coord.offsetBy(1, 0, 0), vel);
        accessor.setValue(coord.offsetBy(1, 0, 1), vel);
        accessor.setValue(coord.offsetBy(1, 1, 0), vel);
        accessor.setValue(coord.offsetBy(1, 1, 1), vel);
    };

    inline void setVelX(Accessor &accessor, int k, int j, int i, double vel) {
        openvdb::Vec3d value = accessor.getValue(openvdb::Coord(k, j, i));
        value[0] = vel;
        accessor.setValue(openvdb::Coord(k, j, i), value);
        accessor.setValue(openvdb::Coord(k, j, i + 1), value);
        accessor.setValue(openvdb::Coord(k, j + 1, i), value);
        accessor.setValue(openvdb::Coord(k, j + 1, i + 1), value);
        accessor.setValue(openvdb::Coord(k + 1, j, i), value);
        accessor.setValue(openvdb::Coord(k + 1, j, i + 1), value);
        accessor.setValue(openvdb::Coord(k + 1, j + 1, i), value);
        accessor.setValue(openvdb::Coord(k + 1, j + 1, i + 1), value);
    };

    inline void setVelX(Accessor &accessor, openvdb::Coord coord, double vel) {
        openvdb::Vec3d value = accessor.getValue(coord);
        value[0] = vel;
        accessor.setValue(coord, value);
        accessor.setValue(coord.offsetBy(0, 0, 1), value);
        accessor.setValue(coord.offsetBy(0, 1, 0), value);
        accessor.setValue(coord.offsetBy(0, 1, 1), value);
        accessor.setValue(coord.offsetBy(1, 0, 0), value);
        accessor.setValue(coord.offsetBy(1, 0, 1), value);
        accessor.setValue(coord.offsetBy(1, 1, 0), value);
        accessor.setValue(coord.offsetBy(1, 1, 1), value);
    };

    inline void setVelY(Accessor &accessor, int k, int j, int i, double vel) {
        openvdb::Vec3d value = accessor.getValue(openvdb::Coord(k, j, i));
        value[1] = vel;
        accessor.setValue(openvdb::Coord(k, j, i), value);
        accessor.setValue(openvdb::Coord(k, j, i + 1), value);
        accessor.setValue(openvdb::Coord(k, j + 1, i), value);
        accessor.setValue(openvdb::Coord(k, j + 1, i + 1), value);
        accessor.setValue(openvdb::Coord(k + 1, j, i), value);
        accessor.setValue(openvdb::Coord(k + 1, j, i + 1), value);
        accessor.setValue(openvdb::Coord(k + 1, j + 1, i), value);
        accessor.setValue(openvdb::Coord(k + 1, j + 1, i + 1), value);
    };

    inline void setVelY(Accessor &accessor, openvdb::Coord coord, double vel) {
        openvdb::Vec3d value = accessor.getValue(coord);
        value[1] = vel;
        accessor.setValue(coord, value);
        accessor.setValue(coord.offsetBy(0, 0, 1), value);
        accessor.setValue(coord.offsetBy(0, 1, 0), value);
        accessor.setValue(coord.offsetBy(0, 1, 1), value);
        accessor.setValue(coord.offsetBy(1, 0, 0), value);
        accessor.setValue(coord.offsetBy(1, 0, 1), value);
        accessor.setValue(coord.offsetBy(1, 1, 0), value);
        accessor.setValue(coord.offsetBy(1, 1, 1), value);
    };

    inline void setVelZ(Accessor &accessor, int k, int j, int i, double vel) {
        openvdb::Vec3d value = accessor.getValue(openvdb::Coord(k, j, i));
        value[2] = vel;
        accessor.setValue(openvdb::Coord(k, j, i), value);
        accessor.setValue(openvdb::Coord(k, j, i + 1), value);
        accessor.setValue(openvdb::Coord(k, j + 1, i), value);
        accessor.setValue(openvdb::Coord(k, j + 1, i + 1), value);
        accessor.setValue(openvdb::Coord(k + 1, j, i), value);
        accessor.setValue(openvdb::Coord(k + 1, j, i + 1), value);
        accessor.setValue(openvdb::Coord(k + 1, j + 1, i), value);
        accessor.setValue(openvdb::Coord(k + 1, j + 1, i + 1), value);
    };

    inline void setVelZ(Accessor &accessor, openvdb::Coord coord, double vel) {
        openvdb::Vec3d value = accessor.getValue(coord);
        value[2] = vel;
        accessor.setValue(coord, value);
        accessor.setValue(coord.offsetBy(0, 0, 1), value);
        accessor.setValue(coord.offsetBy(0, 1, 0), value);
        accessor.setValue(coord.offsetBy(0, 1, 1), value);
        accessor.setValue(coord.offsetBy(1, 0, 0), value);
        accessor.setValue(coord.offsetBy(1, 0, 1), value);
        accessor.setValue(coord.offsetBy(1, 1, 0), value);
        accessor.setValue(coord.offsetBy(1, 1, 1), value);
    };

    inline void addToVelInterpolated(Accessor &accessor, float x, float y, float z, openvdb::Vec3d value) {
        // Calculate indices
        auto isPoint = i2w_transform->worldToIndex(openvdb::Vec3d(x, y, z));
        const int i = isPoint[0];
        const int j = isPoint[1];
        const int k = isPoint[2];

        const double i_frac = isPoint[0] - i;
        const double j_frac = isPoint[1] - j;
        const double k_frac = isPoint[2] - k;

        // Spread in z
        openvdb::Vec3d value_0 = (1 - k_frac) * value;
        openvdb::Vec3d value_1 = k_frac * value;

        // Spread in y
        openvdb::Vec3d value_00 = (1 - j_frac) * value_0;
        openvdb::Vec3d value_10 = j_frac * value_0;
        openvdb::Vec3d value_01 = (1 - j_frac) * value_1;
        openvdb::Vec3d value_11 = j_frac * value_1;

        // Spread in x
        openvdb::Vec3d value_000 = (1 - i_frac) * value_00;
        openvdb::Vec3d value_100 = i_frac * value_00;
        openvdb::Vec3d value_010 = (1 - i_frac) * value_10;
        openvdb::Vec3d value_110 = i_frac * value_10;
        openvdb::Vec3d value_001 = (1 - i_frac) * value_01;
        openvdb::Vec3d value_101 = i_frac * value_01;
        openvdb::Vec3d value_011 = (1 - i_frac) * value_11;
        openvdb::Vec3d value_111 = i_frac * value_11;

        auto prev_v_000 = accessor.getValue(openvdb::Coord(i, j, k));
        auto prev_v_100 = accessor.getValue(openvdb::Coord(i + 1, j, k));
        auto prev_v_010 = accessor.getValue(openvdb::Coord(i, j + 1, k));
        auto prev_v_110 = accessor.getValue(openvdb::Coord(i + 1, j + 1, k));
        auto prev_v_001 = accessor.getValue(openvdb::Coord(i, j, k + 1));
        auto prev_v_101 = accessor.getValue(openvdb::Coord(i + 1, j, k + 1));
        auto prev_v_011 = accessor.getValue(openvdb::Coord(i, j + 1, k + 1));
        auto prev_v_111 = accessor.getValue(openvdb::Coord(i + 1, j + 1, k + 1));


        // Write data
        accessor.setValue(openvdb::Coord(i, j, k), prev_v_000 + value_000);
        accessor.setValue(openvdb::Coord(i + 1, j, k), prev_v_100 + value_100);
        accessor.setValue(openvdb::Coord(i, j + 1, k), prev_v_010 + value_010);
        accessor.setValue(openvdb::Coord(i + 1, j + 1, k), prev_v_110 + value_110);
        accessor.setValue(openvdb::Coord(i, j, k + 1), prev_v_001 + value_001);
        accessor.setValue(openvdb::Coord(i + 1, j, k + 1), prev_v_101 + value_101);
        accessor.setValue(openvdb::Coord(i, j + 1, k + 1), prev_v_011 + value_011);
        accessor.setValue(openvdb::Coord(i + 1, j + 1, k + 1), prev_v_111 + value_111);
    };

    inline void setVelInterpolated(Accessor &accessor, float x, float y, float z, openvdb::Vec3d value) {
        // Calculate indices
        auto isPoint = i2w_transform->worldToIndex(openvdb::Vec3d(x, y, z));
        const int i = isPoint[0];
        const int j = isPoint[1];
        const int k = isPoint[2];

        const double i_frac = isPoint[0] - i;
        const double j_frac = isPoint[1] - j;
        const double k_frac = isPoint[2] - k;

        // Spread in z
        openvdb::Vec3d value_0 = (1 - k_frac) * value;
        openvdb::Vec3d value_1 = k_frac * value;

        // Spread in y
        openvdb::Vec3d value_00 = (1 - j_frac) * value_0;
        openvdb::Vec3d value_10 = j_frac * value_0;
        openvdb::Vec3d value_01 = (1 - j_frac) * value_1;
        openvdb::Vec3d value_11 = j_frac * value_1;

        // Spread in x
        openvdb::Vec3d value_000 = (1 - i_frac) * value_00;
        openvdb::Vec3d value_100 = i_frac * value_00;
        openvdb::Vec3d value_010 = (1 - i_frac) * value_10;
        openvdb::Vec3d value_110 = i_frac * value_10;
        openvdb::Vec3d value_001 = (1 - i_frac) * value_01;
        openvdb::Vec3d value_101 = i_frac * value_01;
        openvdb::Vec3d value_011 = (1 - i_frac) * value_11;
        openvdb::Vec3d value_111 = i_frac * value_11;

        // Write data
        accessor.setValue(openvdb::Coord(i, j, k), value_000);
        accessor.setValue(openvdb::Coord(i + 1, j, k), value_100);
        accessor.setValue(openvdb::Coord(i, j + 1, k), value_010);
        accessor.setValue(openvdb::Coord(i + 1, j + 1, k), value_110);
        accessor.setValue(openvdb::Coord(i, j, k + 1), value_001);
        accessor.setValue(openvdb::Coord(i + 1, j, k + 1), value_101);
        accessor.setValue(openvdb::Coord(i, j + 1, k + 1), value_011);
        accessor.setValue(openvdb::Coord(i + 1, j + 1, k + 1), value_111);
    };

    void swapVelocityBuffers();
    void swapValidMasks();

private:
    openvdb::Vec3dGrid::Ptr _vel_front;
    openvdb::Vec3dGrid::Ptr _vel_back;
    openvdb::Vec3dGrid::Ptr _vel_prev;
    openvdb::Vec3dGrid::Ptr _vel_diff;
    openvdb::math::Transform::Ptr i2w_transform;
    openvdb::MaskGrid::Ptr valid_mask_u_front_buffer;
    openvdb::MaskGrid::Ptr valid_mask_u_back_buffer;
    openvdb::MaskGrid::Ptr valid_mask_v_front_buffer;
    openvdb::MaskGrid::Ptr valid_mask_v_back_buffer;
    openvdb::MaskGrid::Ptr valid_mask_w_front_buffer;
    openvdb::MaskGrid::Ptr valid_mask_w_back_buffer;

    openvdb::FloatGrid::Ptr u_weights;
    openvdb::FloatGrid::Ptr v_weights;
    openvdb::FloatGrid::Ptr w_weights;
};

#endif
