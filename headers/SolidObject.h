#ifndef SOLIDOBJECT_H
#define SOLIDOBJECT_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <LevelSet.h>
#include <openvdb/openvdb.h>

class SolidObject {
public:
    SolidObject(float radius, float voxel_size, double density, openvdb::math::Transform::Ptr xform,
                openvdb::Vec3d center, openvdb::Vec3d vel_translational, openvdb::Vec3d vel_angular);

    inline LevelSet &levelSet() { return _level_set; }
    inline openvdb::FloatGrid::Ptr &uWeights() { return _u_weights; }
    inline openvdb::FloatGrid::Ptr &vWeights() { return _v_weights; }
    inline openvdb::FloatGrid::Ptr &wWeights() { return _w_weights; }
    inline openvdb::Vec3d &velTranslational() { return _vel_translational; }
    inline openvdb::Vec3d &velAngular() { return _vel_angular; }

    inline openvdb::Vec3d centerOfMass() { return _center_of_mass; }
    inline Eigen::SparseMatrix<double, 1> &massMatrixInverse() { return _mass_matrix_inverse; }
    inline Eigen::SparseMatrix<double, 1> &jMatrix() { return _j_matrix; }

    void update_positions_orientations(float dt);
    void compute_face_fractions();
    Eigen::Vector<double, 6> combinedVelocityVector();
    void setVelocityVector(Eigen::Vector< double, 6> vel_vector);
    void generateVelocityField();

    // Constructs the sparse JMatrix
    void constructJMatrix(int idx, openvdb::Coord coord, float voxelSize);

    void add_forces(float dt);

private:
    openvdb::Vec3d _center_of_mass;
    LevelSet _level_set;
    Eigen::SparseMatrix<double, 1> _mass_matrix_inverse;
    openvdb::Vec3d _vel_translational, _vel_angular;
    openvdb::Vec3dGrid::Ptr vel_field;
    openvdb::math::Transform::Ptr xform;
    
    // Density in kg/m^3
    double _density;

    // Solid face area fractions
    openvdb::FloatGrid::Ptr _u_weights = openvdb::createGrid<openvdb::FloatGrid>(0);
    openvdb::FloatGrid::Ptr _v_weights = openvdb::createGrid<openvdb::FloatGrid>(0);
    openvdb::FloatGrid::Ptr _w_weights = openvdb::createGrid<openvdb::FloatGrid>(0);

    // These accessors are used in methods that are called inside iterations
    // (like constructJMatrix in pressure solve) to provide access speedups
    openvdb::FloatGrid::Accessor _u_weight_accessor;
    openvdb::FloatGrid::Accessor _v_weight_accessor;
    openvdb::FloatGrid::Accessor _w_weight_accessor;


    // Structures used for pressure solve
    Eigen::SparseMatrix<double, 1> _j_matrix;
    Eigen::Vector<double, 6> _solid_vel_vector;
};

#endif