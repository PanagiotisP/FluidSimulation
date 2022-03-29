#ifndef LEVELSET_H
#define LEVELSET_H

#include "Grid.h"
#include "util.h"
#include "vec.h"

#include <functional>
#include <iostream>
#include <math.h>

typedef std::function<float(float, float)> DistanceFunction;

class LevelSet: public Grid<float> {
public:
    LevelSet(int size_x, int size_y, float length_x, float length_y);
    ~LevelSet();

    float distance(int from_i, int from_j, int to_i, int to_j);

    float computeUpwindDiffX(int i, int j);
    float computeUpwindDiffY(int i, int j);
    float computeUpwindGradientX(int i, int j, float vel_x);
    float computeUpwindGradientY(int i, int j, float vel_y);
    void construct_from_edges(const std::vector<Vec2f> &points, bool inverse=false);
    void construct_from_points(const std::vector<Vec2f> &points);
    void construct_from_known_geometry(DistanceFunction f);

    void union_level_set(LevelSet & level_set);
    void intersect_level_set(LevelSet & level_set);

    void add_to_set(float val);
    // Takes an inside-outside level set that is far from sdf and produces sdf
    void redistance();

    LevelSet invert();
    friend std::ostream &operator<<(std::ostream &out, LevelSet &l);

private:
    // Eikonal PDE for distance solver
    // Propagates a level set that is defined correctly near the surface boundary to the whole grid
    // Markers indicate which grid cells have fixed values from initialisation
    // Cross markers indicate which grid edges have crossings, that help us to correctly evaluate the inside-outside
    void fast_sweeping(Grid<int> &markers, bool isSigned, int max_iters = 2);
    void transform_to_signed_distance(Grid<int> &cross_markers, bool inverse);
};

#endif