#pragma once
#ifndef RENDER_UTILS_H
#define RENDER_UTILS_H

#include "FluidSimulator.h"
#include "LevelSet.h"
#include "MacGrid.h"

#include <SFML/Graphics.hpp>
#include <array>
#include <vector>

// Axis aligned rectangular box
template <typename T>
DistanceFunction distance_from_axis_aligned_box(BBox<T> box_coords, bool inverted = false) {
    return [&box_coords, inverted](T x, T y) {
        T distance;
        if (box_coords.x_min < x && x < box_coords.x_max && box_coords.y_min < y && y < box_coords.y_max) {
            distance = max(box_coords.x_min - x, x - box_coords.x_max, box_coords.y_min - y, y - box_coords.y_max);
        } else {
            T closest_x, closest_y;
            if (x < box_coords.x_min)
                closest_x = box_coords.x_min;
            else if (box_coords.x_max < x)
                closest_x = box_coords.x_max;
            else
                closest_x = x;
            if (y < box_coords.y_min)
                closest_y = box_coords.y_min;
            else if (box_coords.y_max < y)
                closest_y = box_coords.y_max;
            else
                closest_y = y;
            distance = sqrt(sqr(closest_x - x) + sqr(closest_y - y));
        }
        return inverted ? -distance : distance;
    };
}

// Sphere
template <typename T>
DistanceFunction distance_from_sphere(T center_x, T center_y, T radius, bool inverted = false) {
    return [center_x, center_y, radius, inverted](T x, T y) {
        T distance = dist(Vec<2, T>(center_x, center_y), Vec<2, T>(x, y)) - radius;
        return inverted ? -distance : distance;
    };
}

#endif