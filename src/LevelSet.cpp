#include "LevelSet.h"

LevelSet::LevelSet(int size_x, int size_y, float length_x, float length_y):
 Grid<float>(size_x, size_y, length_x / size_x, length_y / size_y) {
    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) { (*this)(i, j) = 0; }
    }
}

LevelSet::~LevelSet() {}

float LevelSet::distance(int from_i, int from_j, int to_i, int to_j) {
    float i_diff = (to_i - from_i);
    float j_diff = (to_j - from_j);
    return sqrt(i_diff * i_diff + j_diff * j_diff);
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