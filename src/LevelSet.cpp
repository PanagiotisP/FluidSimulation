#include "LevelSet.h"

#include "util.h"

#include <Eigen/Geometry>
#include <functional>
#include <iomanip>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/ParticlesToLevelSet.h>

LevelSet::LevelSet(openvdb::FloatGrid::Ptr level_set): _level_set(level_set) {}

LevelSet::~LevelSet() {}

float LevelSet::valueInterpolatedW(LevelSet::BoxSampler &sampler, float x, float y, float z) {
    // Need to check whether we are thread-safe or not
    return sampler.wsSample(openvdb::Vec3R(x, y, z));
}
float LevelSet::valueInterpolatedW(LevelSet::BoxSampler &sampler, openvdb::Vec3R wsPoint) {
    // Need to check whether we are thread-safe or not
    return sampler.wsSample(wsPoint);
}
float LevelSet::valueInterpolatedI(LevelSet::BoxSampler &sampler, float x, float y, float z) {
    // Need to check whether we are thread-safe or not
    return sampler.isSample(openvdb::Vec3R(x, y, z));
}
float LevelSet::valueInterpolatedI(LevelSet::BoxSampler &sampler, openvdb::Vec3R isPoint) {
    // Need to check whether we are thread-safe or not
    return sampler.isSample(isPoint);
}

void LevelSet::construct_from_points(ParticleSet &points, float voxel_size) {
    _level_set->clear();
    openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*_level_set);
    raster.setRmin(points.getRadius());
    raster.rasterizeSpheres(points, points.getRadius());
    raster.finalize(true);
}

void LevelSet::unionLevelSet(LevelSet &level_set) {
    _level_set = openvdb::tools::csgUnionCopy(*_level_set, *level_set.getLevelSet());
}

void LevelSet::intersectionLevelSet(LevelSet &level_set) {
    _level_set = openvdb::tools::csgIntersectionCopy(*_level_set, *level_set.getLevelSet());
}

void LevelSet::differenceLevelSet(LevelSet &level_set) {
    _level_set = openvdb::tools::csgDifferenceCopy(*_level_set, *level_set.getLevelSet());
}

void LevelSet::cycle_array(double *arr, int size) {
    double t = arr[0];
    for (int i = 0; i < size - 1; ++i) arr[i] = arr[i + 1];
    arr[size - 1] = t;
}

// Given two signed distance values (line endpoints), determine what fraction of a connecting segment is "inside"
double LevelSet::fraction_inside(double phi_left, double phi_right) {
    if (phi_left < 0 && phi_right < 0) return 1;
    if (phi_left < 0 && phi_right >= 0) return phi_left / (phi_left - phi_right);
    if (phi_left >= 0 && phi_right < 0)
        return phi_right / (phi_right - phi_left);
    else
        return 0;
}


// Given four signed distance values (square corners), determine what fraction of the square is "inside"
double LevelSet::fraction_inside(double phi_bl, double phi_br, double phi_tl, double phi_tr) {

    int inside_count = (phi_bl < 0 ? 1 : 0) + (phi_tl < 0 ? 1 : 0) + (phi_br < 0 ? 1 : 0) + (phi_tr < 0 ? 1 : 0);
    double list[] = { phi_bl, phi_br, phi_tr, phi_tl };
    double return_value;
    if (inside_count == 4)
        return 1;
    else if (inside_count == 3) {
        // rotate until the positive value is in the first position
        while (list[0] < 0) { cycle_array(list, 4); }

        // Work out the area of the exterior triangle
        double side0 = 1 - fraction_inside(list[0], list[3]);
        double side1 = 1 - fraction_inside(list[0], list[1]);
        return_value = 1 - 0.5f * side0 * side1;
    } else if (inside_count == 2) {

        // rotate until a negative value is in the first position, and the next negative is in either slot 1 or 2.
        while (list[0] >= 0 || !(list[1] < 0 || list[2] < 0)) { cycle_array(list, 4); }

        if (list[1] < 0) { // the matching signs are adjacent
            double side_left = fraction_inside(list[0], list[3]);
            double side_right = fraction_inside(list[1], list[2]);
            return_value = 0.5f * (side_left + side_right);
        } else { // matching signs are diagonally opposite
            // determine the centre point's sign to disambiguate this case
            double middle_point = 0.25f * (list[0] + list[1] + list[2] + list[3]);
            if (middle_point < 0) {
                double area = 0;

                // first triangle (top left)
                double side1 = 1 - fraction_inside(list[0], list[3]);
                double side3 = 1 - fraction_inside(list[2], list[3]);

                area += 0.5f * side1 * side3;

                // second triangle (top right)
                double side2 = 1 - fraction_inside(list[2], list[1]);
                double side0 = 1 - fraction_inside(list[0], list[1]);
                area += 0.5f * side0 * side2;

                return_value = 1 - area;
            } else {
                double area = 0;

                // first triangle (bottom left)
                double side0 = fraction_inside(list[0], list[1]);
                double side1 = fraction_inside(list[0], list[3]);
                area += 0.5f * side0 * side1;

                // second triangle (top right)
                double side2 = fraction_inside(list[2], list[1]);
                double side3 = fraction_inside(list[2], list[3]);
                area += 0.5f * side2 * side3;
                return_value = area;
            }
        }
    } else if (inside_count == 1) {
        // rotate until the negative value is in the first position
        while (list[0] >= 0) { cycle_array(list, 4); }

        // Work out the area of the interior triangle, and subtract from 1.
        double side0 = fraction_inside(list[0], list[3]);
        double side1 = fraction_inside(list[0], list[1]);
        return_value = 0.5f * side0 * side1;
    } else
        return_value = 0;

    if (return_value > 0.9)
        return 1;
    else
        return clamp(return_value, 0., 1.);
}