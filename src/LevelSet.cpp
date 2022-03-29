#include "LevelSet.h"

#include "util.h"

#include <Eigen/Geometry>
#include <functional>
#include <iomanip>
LevelSet::LevelSet(int size_x, int size_y, float length_x, float length_y):
 Grid<float>(size_x, size_y, length_x / size_x, length_y / size_y) {
    float max_val = size_x * size_y;
    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) { (*this)(i, j) = max_val; }
    }
}

LevelSet::~LevelSet() {}

float LevelSet::distance(int from_i, int from_j, int to_i, int to_j) {
    float i_diff = (to_i - from_i);
    float j_diff = (to_j - from_j);
    return sqrt(i_diff * i_diff + j_diff * j_diff);
}

float LevelSet::computeUpwindDiffX(int i, int j) {
    int i_plus_1 = i + 1;
    int i_minus_1 = i - 1;
    // Border cases
    i_plus_1 = clamp(i_plus_1, 0, _SIZE_X - 1);
    i_minus_1 = clamp(i_minus_1, 0, _SIZE_X - 1);

    if (i_minus_1 < 0)
        return (*this)(i_plus_1, j);
    else if (i_plus_1 > _SIZE_X - 1)
        return (*this)(i_minus_1, j);
    else
        return min((*this)(i_minus_1, j), (*this)(i_plus_1, j));
}

float LevelSet::computeUpwindDiffY(int i, int j) {
    int j_plus_1 = j + 1;
    int j_minus_1 = j - 1;
    // Border cases
    j_plus_1 = clamp(j_plus_1, 0, _SIZE_Y - 1);
    j_minus_1 = clamp(j_minus_1, 0, _SIZE_Y - 1);

    if (j_minus_1 < 0)
        return (*this)(i, j_plus_1);
    else if (j_plus_1 > _SIZE_X - 1)
        return (*this)(i, j_minus_1);
    else
        return min((*this)(i, j_minus_1), (*this)(i, j_plus_1));
}

float LevelSet::computeUpwindGradientX(int i, int j, float vel_x) {
    int i_plus1 = i + 1;
    int i_minus1 = i - 1;
    // Border cases
    i = clamp(i, 0, _SIZE_X - 1);
    j = clamp(j, 0, _SIZE_Y - 1);
    i_plus1 = clamp(i_plus1, 0, _SIZE_X - 1);
    i_minus1 = clamp(i_minus1, 0, _SIZE_X - 1);

    float grad_x = vel_x < 0 ? (*this)(i_plus1, j) - (*this)(i, j) : (*this)(i, j) - (*this)(i_minus1, j);
    return grad_x;
}

float LevelSet::computeUpwindGradientY(int i, int j, float vel_y) {
    int j_plus1 = j + 1;
    int j_minus1 = j - 1;
    // Border cases
    i = clamp(i, 0, _SIZE_X - 1);
    j = clamp(j, 0, _SIZE_Y - 1);
    j_plus1 = clamp(j_plus1, 0, _SIZE_Y - 1);
    j_minus1 = clamp(j_minus1, 0, _SIZE_Y - 1);

    float grad_y = vel_y < 0 ? (*this)(i, j_plus1) - (*this)(i, j) : (*this)(i, j) - (*this)(i, j_minus1);
    return grad_y;
}

// Initialise level set from a set of points that indicate edges.
// Subsequent points are assumed to be connected by straight lines.
// The last point is connected to the first one to produce a watertight surface.
void LevelSet::construct_from_edges(const std::vector<Vec2f> &points, bool inverse) {
    LevelSet &l = *this;
    Grid<int> markers(_SIZE_X, _SIZE_Y, _DELTA_X, _DELTA_Y);
    Grid<int> cross_markers(_SIZE_X, _SIZE_Y, _DELTA_X, _DELTA_Y);

    // Initialise markers to -1 and crosses to 0
    for (int j = 0; j < markers.sizeY(); ++j) {
        for (int i = 0; i < markers.sizeX(); ++i) {
            markers(i, j) = -1;
            cross_markers(i, j) = 0;
        }
    }
    // Set markers and distances for fixed points
    for (int k = 0; k < points.size(); ++k) {
        // Last point is connected with the first one
        int k_next = (k + 1) % points.size();
        Eigen::Vector2<float> start((points[k][0] - 0.5 * _DELTA_X) / _DELTA_X,
                                    (points[k][1] - 0.5 * _DELTA_Y) / _DELTA_Y);
        Eigen::Vector2<float> end((points[k_next][0] - 0.5 * _DELTA_X) / _DELTA_X,
                                  (points[k_next][1] - 0.5 * _DELTA_Y) / _DELTA_Y);
        auto edge(Eigen::Hyperplane<float, 2>::Through(start, end));

        // Implementation of modified Bresehnam's line algorithm from https://gamedev.stackexchange.com/a/182143
        int i = floor(start[0]);
        int j = floor(start[1]);

        float diffX = end[0] - start[0];
        float diffY = end[1] - start[1];

        // Assumes that grid spacing is 1
        float stepX = diffX > 0 ? 1 : -1;
        float stepY = diffY > 0 ? 1 : -1;


        float xOffset = end[0] > start[0] ? (ceil(start[0]) - start[0]) : (start[0] - floor(start[0]));
        float yOffset = end[1] > start[1] ? (ceil(start[1]) - start[1]) : (start[1] - floor(start[1]));
        // Instead of using cos, sin, atan, we will use the following formulas
        // cos(atan(x)) = 1 / sqrt(x^2 + 1)
        // sin(atan(x)) = x / sqrt(x^2 + 1)
        // float angle = atan2(-diffY, diffX);
        // float tMaxX = xOffset / cos(angle);
        // float tMaxY = yOffset / sin(angle);
        // float tDeltaX = 1.0 / cos(angle);
        // float tDeltaY = 1.0 / sin(angle);

        // If denominator is zero set tan to something big, instead of dealing with exceptions
        float tan;
        if (diffX == 0)
            tan = 1e15;
        else if (diffY == 0)
            tan = 1e-15;
        else
            tan = diffY / diffX;
        // How far to move along the ray to cross the first vertical grid cell boundary.
        float tMaxX = diffX != 0 ? xOffset * sqrt(sqr(tan) + 1) : 1e15;
        // How far to move along the ray to cross the first horizontal grid cell boundary.
        float tMaxY = diffY != 0 ? yOffset * sqrt(sqr(tan) + 1) / tan : 1e15;
        // How far to move along the ray to move horizontally 1 grid cell.
        float tDeltaX = 1.0 * sqrt(sqr(tan) + 1);
        // How far to move along the ray to move vertically 1 grid cell.
        float tDeltaY = 1.0 * sqrt(sqr(tan) + 1) / tan;

        // Travel one grid cell at a time.
        float manhattanDistance = fabs(floor(end[0]) - floor(start[0])) + fabs(floor(end[1]) - floor(start[1]));
        for (int t = 0; t <= manhattanDistance; ++t) {
            // If ray is inside a cell of the grid calculate distance for that cell and four adjacent ones
            if (0 <= i && i < _SIZE_X && 0 <= j && j < _SIZE_Y) {
                // Four adjacent cells in cross pattern (e.g. for (4,4) calculate distance for (3,4), (5,4), (4,3),
                // (4,5)
                for (int o = -1; o <= 1; o += 2) {
                    int i_valid = clamp(i + o, 0, _SIZE_X - 1);
                    int j_valid = clamp(j + o, 0, _SIZE_X - 1);

                    markers(i_valid, j) = 1;
                    markers(i, j_valid) = 1;

                    l(i, j_valid) = min(l(i, j_valid), edge.absDistance(Eigen::Vector2<float>(i, j_valid)));
                    l(i_valid, j) = min(l(i_valid, j), edge.absDistance(Eigen::Vector2<float>(i_valid, j)));
                }
                markers(i, j) = 1;
                l(i, j) = min(l(i, j), edge.absDistance(Eigen::Vector2<float>(i, j)));
            }
            // Only move in either X or Y coordinates, not both.
            if (fabs(tMaxX) < fabs(tMaxY)) {
                tMaxX += tDeltaX;
                i += stepX;
            } else {
                if (0 <= i && i < _SIZE_X && 0 <= j && j < _SIZE_Y) cross_markers(i, j)++;
                tMaxY += tDeltaY;
                j += stepY;
            }
        }
    }
    fast_sweeping(markers, false);
    transform_to_signed_distance(cross_markers, inverse);
}

void LevelSet::construct_from_points(const std::vector<Vec2f> &points) {
    LevelSet l(_SIZE_X, _SIZE_Y, _SIZE_X * _DELTA_X, _SIZE_Y * _DELTA_Y);
    Grid<int> markers(_SIZE_X, _SIZE_Y, _DELTA_X, _DELTA_Y);

    // Initialise markers to -1
    for (int j = 0; j < markers.sizeY(); ++j) {
        for (int i = 0; i < markers.sizeX(); ++i) { markers(i, j) = -1; }
    }
    // Set markers and distances for fixed points
    for (auto point : points) {
        Vec2f point_scaled((point[0] - 0.5 * _DELTA_X) / _DELTA_X, (point[1] - 0.5 * _DELTA_Y) / _DELTA_Y);
        int i_min = clamp(static_cast<int>(floor(point_scaled[0])), 0, _SIZE_X - 1);
        int j_min = clamp(static_cast<int>(floor(point_scaled[1])), 0, _SIZE_Y - 1);
        int i_max = clamp(static_cast<int>(ceil(point_scaled[0])), 0, _SIZE_X - 1);
        int j_max = clamp(static_cast<int>(ceil(point_scaled[1])), 0, _SIZE_Y - 1);
        for (int j = j_min; j <= j_max; ++j) {
            for (int i = i_min; i <= i_max; ++i) {
                markers(i, j) = 1;
                l(i, j) = min(l(i, j), dist(point_scaled, Vec2f(i, j)));
            }
        }
    }
    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) {
            if (markers(i, j) == -1) l(i, j) = INFINITY;
        }
    }
    *this = l;
    fast_sweeping(markers, false);
}

void LevelSet::construct_from_known_geometry(DistanceFunction f) {
    auto &l = *this;
    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) { 
            l(i, j) = f((i + 0.f) * _DELTA_X, (j + 0.f) * _DELTA_Y); }
    }
}

void LevelSet::fast_sweeping(Grid<int> &markers, bool isSigned, int max_iters) {
    LevelSet &l = *this;
    int grid_size = 1;
    Grid<float> sign(_SIZE_X, _SIZE_Y, _DELTA_X, _DELTA_Y);
    if (isSigned) {
        for (int j = 0; j < _SIZE_Y; ++j) {
            for (int i = 0; i < _SIZE_X; ++i) {
                sign(i, j) = sgn(l(i, j));
                l(i, j) = fabs(l(i, j));
            }
        }
    }

    std::function<void(int, int)> sweepOne = [&l, grid_size](int i, int j) {
        float phi_0 = l.computeUpwindDiffX(i, j);
        float phi_1 = l.computeUpwindDiffY(i, j);
        sort(phi_0, phi_1);

        float unique_solution = phi_0 + grid_size;

        if (unique_solution > phi_1) {
            unique_solution = 0.5 * (phi_0 + phi_1 + sqrt(2 * grid_size * grid_size - sqr(phi_1 - phi_0)));
        }
        l(i, j) = min(l(i, j), unique_solution);
    };


    for (int iter = 0; iter < max_iters; ++iter) {
        for (int j = 0; j < _SIZE_Y; ++j) {
            for (int i = 0; i < _SIZE_X; ++i) {
                if (markers(i, j) == -1) { 
                    sweepOne(i, j); }
            }
        }

        for (int j = 0; j < _SIZE_Y; ++j) {
            for (int i = _SIZE_X - 1; i >= 0; --i) {
                if (markers(i, j) == -1) { sweepOne(i, j); }
            }
        }

        for (int j = _SIZE_Y - 1; j >= 0; --j) {
            for (int i = 0; i < _SIZE_X; ++i) {
                if (markers(i, j) == -1) { sweepOne(i, j); }
            }
        }

        for (int j = _SIZE_Y - 1; j >= 0; --j) {
            for (int i = _SIZE_X - 1; i >= 0; --i) {
                if (markers(i, j) == -1) { sweepOne(i, j); }
            }
        }
    }

    if (isSigned) {
        for (int j = 0; j < _SIZE_Y; ++j) {
            for (int i = 0; i < _SIZE_X; ++i) { l(i, j) *= sign(i, j); }
        }
    }
}

void LevelSet::transform_to_signed_distance(Grid<int> &cross_markers, bool inverse) {
    auto &l = *this;
    // Fix sign based on crosses
    for (int j = 0; j < _SIZE_Y; ++j) {
        int intersection_count = 0;
        for (int i = 0; i < _SIZE_X; ++i) {
            if(inverse) l(i,j) *= -1;
            if (intersection_count % 2 == 1) { l(i, j) *= -1; }
            intersection_count += cross_markers(i, j);
        }
    }
}

void LevelSet::union_level_set(LevelSet &level_set) {
    for (int j = 0; j < (*this).sizeY(); ++j) {
        for (int i = 0; i < (*this).sizeX(); ++i) { (*this)(i, j) = min((*this)(i, j), level_set(i, j)); }
    }
}

void LevelSet::intersect_level_set(LevelSet &level_set) {
    for (int j = 0; j < (*this).sizeY(); ++j) {
        for (int i = 0; i < (*this).sizeX(); ++i) { (*this)(i, j) = max((*this)(i, j), level_set(i, j)); }
    }
}

LevelSet LevelSet::invert() {
    LevelSet new_level_set(*this);
    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) { new_level_set(i, j) = -1 * (*this)(i, j); }
    }
    return new_level_set;
}

std::ostream &operator<<(std::ostream &out, LevelSet &l) {
    for (int j = l.sizeY() - 1; j >= 0; --j) {
        out << std::endl;
        for (int i = 0; i < l.sizeX(); ++i) { out << std::setw(5) << l(i, j) << ' '; }
    }
    out << std::endl;
    return out;
}

void LevelSet::add_to_set(float val) {
    val = val / min(_DELTA_X, _DELTA_Y);
    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) { (*this)(i, j) += val; }
    }
}

void LevelSet::redistance() {
    auto &l = (*this);
    LevelSet l_prev(l);
    Grid<int> markers(_SIZE_X, _SIZE_Y, _DELTA_X, _DELTA_Y);
    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) { markers(i, j) = -1; }
    }

    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) {
            int sign = sgn(l_prev(i, j));
            int i_minus_1 = clamp(i - 1, 0, _SIZE_X - 1);
            int j_minus_1 = clamp(j - 1, 0, _SIZE_Y - 1);
            int i_plus_1 = clamp(i + 1, 0, _SIZE_X - 1);
            int j_plus_1 = clamp(j + 1, 0, _SIZE_Y - 1);

            // If at least one neighbour has different sign then the point is close to the surface
            if (sign != sgn(l_prev(i_minus_1, j)) || sign != sgn(l_prev(i, j_minus_1)) || sign != sgn(l_prev(i_plus_1, j))
                || sign != sgn(l_prev(i, j_plus_1))) {
                markers(i, j) = 1;

                float cur_val = l_prev(i, j);
                l(i, j) = fabs(l(i, j));
                if (sign != sgn(l_prev(i_minus_1, j)))
                    l(i, j) = min(l(i, j), fabs(cur_val / (cur_val - l_prev(i_minus_1, j))));
                if (sign != sgn(l_prev(i, j_minus_1)))
                    l(i, j) = min(l(i, j), fabs(cur_val / (cur_val - l_prev(i, j_minus_1))));
                if (sign != sgn(l_prev(i_plus_1, j)))
                    l(i, j) = min(l(i, j), fabs(cur_val / (cur_val - l_prev(i_plus_1, j))));
                if (sign != sgn(l_prev(i, j_plus_1)))
                    l(i, j) = min(l(i, j), fabs(cur_val / (cur_val - l_prev(i, j_plus_1))));
                l(i, j) *= sign;
            }
        }
    }

    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) {
            if (markers(i, j) == -1) l(i, j) = sgn(l(i, j)) * INFINITY;
        }
    }
    fast_sweeping(markers, true);
}