#include "MacGrid.h"

MacGrid::MacGrid(openvdb::math::Transform::Ptr i2w_transform, float voxel_size): _i2w_transform(i2w_transform), _voxel_size(voxel_size) {
    _vel_front = openvdb::createGrid<openvdb::Vec3dGrid>(openvdb::Vec3d(0));
    _vel_front->setTransform(i2w_transform);
    _vel_front->setGridClass(openvdb::GRID_STAGGERED);
    _vel_back = openvdb::createGrid<openvdb::Vec3dGrid>(openvdb::Vec3d(0));
    _vel_back->setTransform(i2w_transform);
    _vel_back->setGridClass(openvdb::GRID_STAGGERED);
    _vel_prev = openvdb::createGrid<openvdb::Vec3dGrid>(openvdb::Vec3d(0));
    _vel_prev->setTransform(i2w_transform);
    _vel_prev->setGridClass(openvdb::GRID_STAGGERED);
    _vel_diff = openvdb::createGrid<openvdb::Vec3dGrid>(openvdb::Vec3d(0));
    _vel_diff->setTransform(i2w_transform);
    _vel_diff->setGridClass(openvdb::GRID_STAGGERED);

    _valid_mask_u_front_buffer = openvdb::MaskGrid::create();
    _valid_mask_u_front_buffer->setTransform(i2w_transform);
    _valid_mask_u_back_buffer = openvdb::MaskGrid::create();
    _valid_mask_u_back_buffer->setTransform(i2w_transform);

    _valid_mask_v_front_buffer = openvdb::MaskGrid::create();
    _valid_mask_v_front_buffer->setTransform(i2w_transform);
    _valid_mask_v_back_buffer = openvdb::MaskGrid::create();
    _valid_mask_v_back_buffer->setTransform(i2w_transform);

    _valid_mask_w_front_buffer = openvdb::MaskGrid::create();
    _valid_mask_w_front_buffer->setTransform(i2w_transform);
    _valid_mask_w_back_buffer = openvdb::MaskGrid::create();
    _valid_mask_w_back_buffer->setTransform(i2w_transform);

    _u_weights = openvdb::createGrid<openvdb::FloatGrid>(0);
    _u_weights->setTransform(i2w_transform);
    _v_weights = openvdb::createGrid<openvdb::FloatGrid>(0);
    _v_weights->setTransform(i2w_transform);
    _w_weights = openvdb::createGrid<openvdb::FloatGrid>(0);
    _w_weights->setTransform(i2w_transform);
}

MacGrid::~MacGrid() {}

void MacGrid::updatePreviousVelocityBuffer() { _vel_prev = _vel_front->deepCopy(); }

struct Local {
    static inline void diff(const openvdb::Vec3d &a, const openvdb::Vec3d &b, openvdb::Vec3d &result) {
        result = a - b;
    }
};

void MacGrid::updateDiffBuffers() { _vel_diff->tree().combine2(_vel_front->tree(), _vel_prev->tree(), Local::diff); }

void MacGrid::swapVelocityBuffers() { _vel_front.swap(_vel_back); }
void MacGrid::swapValidMasks() {
    _valid_mask_u_front_buffer.swap(_valid_mask_u_back_buffer);
    _valid_mask_v_front_buffer.swap(_valid_mask_v_back_buffer);
    _valid_mask_w_front_buffer.swap(_valid_mask_w_back_buffer);
}