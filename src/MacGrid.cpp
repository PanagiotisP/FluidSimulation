#include "MacGrid.h"

MacGrid::MacGrid(openvdb::math::Transform::Ptr i2w_transform): i2w_transform(i2w_transform) {
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