#include "MacGrid.h"
#define SOLID_BOUNDARIES

MacGrid::MacGrid(int size_x, int size_y, float length_x, float length_y, float ambient_temperature,
                 float ambient_concentration):
 GridInterface(size_x, size_y, length_x / size_x, length_y / size_y),
 _vel_x_previous(size_x, size_y, _DELTA_X, _DELTA_Y), _vel_y_previous(size_x, size_y, _DELTA_X, _DELTA_Y),
 temperature_previous(size_x, size_y, _DELTA_X, _DELTA_Y), concentration_previous(size_x, size_y, _DELTA_X, _DELTA_Y),
 _vel_x_diff(size_x, size_y, _DELTA_X, _DELTA_Y), _vel_y_diff(size_x, size_y, _DELTA_X, _DELTA_Y),
 temperature_diff(size_x, size_y, _DELTA_X, _DELTA_Y), concentration_diff(size_x, size_y, _DELTA_X, _DELTA_Y),
 _cell_type_buffer(size_x, size_y) {
    _vel_x_front_buffer = std::make_unique< Grid<float> >(size_x, size_y, _DELTA_X, _DELTA_Y);
    _vel_y_front_buffer = std::make_unique< Grid<float> >(size_x, size_y, _DELTA_X, _DELTA_Y);
    _temperature_front_buffer = std::make_unique< Grid<float> >(size_x, size_y, _DELTA_X, _DELTA_Y);
    _concentration_front_buffer = std::make_unique< Grid<float> >(size_x, size_y, _DELTA_X, _DELTA_Y);
    _vel_x_back_buffer = std::make_unique< Grid<float> >(size_x, size_y, _DELTA_X, _DELTA_Y);
    _vel_y_back_buffer = std::make_unique< Grid<float> >(size_x, size_y, _DELTA_X, _DELTA_Y);
    _temperature_back_buffer = std::make_unique< Grid<float> >(size_x, size_y, _DELTA_X, _DELTA_Y);
    _concentration_back_buffer = std::make_unique< Grid<float> >(size_x, size_y, _DELTA_X, _DELTA_Y);

    // Initialise scalars
    for (int j = 0; j < size_y; ++j) {
        for (int i = 0; i < size_x; ++i) {
            (*_temperature_front_buffer)(i, j) = ambient_temperature;
            (*_concentration_front_buffer)(i, j) = ambient_concentration;
        }
    }

    clearCellTypeBuffer();
}

MacGrid::~MacGrid() {}

void MacGrid::clearCellTypeBuffer() {
#ifdef SOLID_BOUNDARIES
    // Set all to SOLID
    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) { _cell_type_buffer(i, j) = SOLID; }
    }
    // Top AIR
    // for (int i = 1; i < _SIZE_X - 1; ++i) _cell_type_buffer(i, _SIZE_Y - 1) = AIR;
    // Set center to LIQUID
    for (int j = 1; j < _SIZE_Y - 1; ++j) {
        for (int i = 1; i < _SIZE_X - 1; ++i) { _cell_type_buffer(i, j) = LIQUID; }
    }
#endif
#ifndef SOLID_BOUNDARIES
    // Set all to AIR
    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) { _cell_type_buffer(i, j) = AIR; }
    }
    // Set center to LIQUID
    for (int j = 1; j < _SIZE_Y - 1; ++j) {
        for (int i = 1; i < _SIZE_X - 1; ++i) { _cell_type_buffer(i, j) = LIQUID; }
    }
#endif // !SOLID_BOUNDARIES
}

void MacGrid::updatePreviousVelocityBuffer() {
    _vel_x_previous = *_vel_x_front_buffer;
    _vel_y_previous = *_vel_y_front_buffer;
    temperature_previous = *_temperature_front_buffer;
    concentration_previous = *_concentration_front_buffer;
}

void MacGrid::updateDiffBuffers() {
    for (int j = 0; j < _SIZE_Y; ++j) {
        for (int i = 0; i < _SIZE_X; ++i) {
            _vel_x_diff(i, j) = _vel_x_front_buffer->value(i, j) - _vel_x_previous.value(i, j);
            _vel_y_diff(i, j) = _vel_y_front_buffer->value(i, j) - _vel_y_previous.value(i, j);
            temperature_diff(i, j) = _temperature_front_buffer->value(i, j) - temperature_previous.value(i, j);
            concentration_diff(i, j) = _concentration_front_buffer->value(i, j) - concentration_previous.value(i, j);
        }
    }
}

// Getters
/*
glm::dmat2 MacGrid::computeVelocityGradientMatrix(int i, int j)
{
    glm::dmat2 vel_grad;
    vel_grad[0][0] = divVelX(i, j);
    vel_grad[0][1] = (_vel_x_front_buffer->value(i, j + 1)
        - _vel_x_front_buffer->value(i, j))
        / _DELTA_Y;
    vel_grad[1][0] = (_vel_y_front_buffer->value(i + 1, j)
        - _vel_y_front_buffer->value(i, j))
        / _DELTA_X;
    vel_grad[1][1] = divVelY(i, j);
    return vel_grad;
}
*/

void MacGrid::swapVelocityBuffers() {
    _vel_x_front_buffer.swap(_vel_x_back_buffer);
    _vel_y_front_buffer.swap(_vel_y_back_buffer);
}


void MacGrid::swapTemperatureBuffers() { _temperature_front_buffer.swap(_temperature_back_buffer); }

void MacGrid::swapConcentrationBuffers() { _concentration_front_buffer.swap(_concentration_back_buffer); }