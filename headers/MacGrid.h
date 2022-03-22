#pragma once

#ifndef MAC_GRID_H
    #define MAC_GRID_H

    #include "Grid.h"

//#include <glm/glm.hpp>

    #include "util.h"

    #include <assert.h>
    #include <iostream>
    #include <memory>

enum CellType { LIQUID, AIR, SOLID };

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


class MacGrid: public GridInterface {
    public:
    MacGrid(int size_x, int size_y, float length_x, float length_y, float ambient_temperature,
            float ambient_concentration);
    ~MacGrid();

    void clearCellTypeBuffer();
    void updatePreviousVelocityBuffer();
    void updateDiffBuffers();

    // Getters
    inline float temperature(int i, int j) const { return _temperature_front_buffer->value(i, j); };

    inline float temperatureBackBuffer(int i, int j) const { return _temperature_back_buffer->value(i, j); };

    inline float concentration(int i, int j) const { return _concentration_front_buffer->value(i, j); };

    inline float concentrationBackBuffer(int i, int j) const { return _concentration_back_buffer->value(i, j); };

    inline float velX(int i, int j) const {
        return (_vel_x_front_buffer->value(i, j) + _vel_x_front_buffer->value(i + 1, j)) / 2;
    };

    inline float velY(int i, int j) const {
        return (_vel_y_front_buffer->value(i, j) + _vel_y_front_buffer->value(i, j + 1)) / 2;
    };

    inline float velXHalfIndexed(int i, int j) const { return _vel_x_front_buffer->value(i, j); };

    inline float velYHalfIndexed(int i, int j) const { return _vel_y_front_buffer->value(i, j); };

    inline float velXBackBufferHalfIndexed(int i, int j) const { return _vel_x_back_buffer->value(i, j); };

    inline float velYBackBufferHalfIndexed(int i, int j) const { return _vel_y_back_buffer->value(i, j); };

    inline float velXBackBuffer(int i, int j) const {
        return (_vel_x_back_buffer->value(i, j) + _vel_x_back_buffer->value(i + 1, j)) / 2;
    };

    inline float velYBackBuffer(int i, int j) const {
        return (_vel_y_back_buffer->value(i, j) + _vel_y_back_buffer->value(i, j + 1)) / 2;
    };

    inline float velXInterpolated(float x, float y) const {
        float v_x = _vel_x_front_buffer->valueInterpolated(x, y - _DELTA_Y * 0.5);
        // -0.5 Due to the MAC grid structure
        return v_x;
    };

    inline float velYInterpolated(float x, float y) const {
        float v_y = _vel_y_front_buffer->valueInterpolated(x - _DELTA_X * 0.5, y);
        // -0.5 Due to the MAC grid structure
        return v_y;
    };

    inline float temperatureInterpolated(float x, float y) const {
        // The subtraction on x and y is necessary due to valueInterpolated function working for the staggered grid
        float temperature = _temperature_front_buffer->valueInterpolated(x - 0.5 * _DELTA_X, y - 0.5 * _DELTA_Y);
        return temperature;
    };

    inline float concentrationInterpolated(float x, float y) const {
        // The subtraction on x and y is necessary due to valueInterpolated function working for the staggered grid
        float concentration = _concentration_front_buffer->valueInterpolated(x - 0.5 * _DELTA_X, y - 0.5 * _DELTA_Y);
        return concentration;
    };

    inline float velXDiffInterpolated(float x, float y) const {
        float v_x = _vel_x_diff.valueInterpolated(x, y - _DELTA_Y * 0.5);
        // -0.5 Due to the MAC grid structure
        return v_x;
    };

    inline float velYDiffInterpolated(float x, float y) const {
        float v_y = _vel_y_diff.valueInterpolated(x - _DELTA_X * 0.5, y);
        // -0.5 Due to the MAC grid structure
        return v_y;
    };

    inline float temperatureDiffInterpolated(float x, float y) const {
        float temperature = temperature_diff.valueInterpolated(x - _DELTA_X * 0.5, y - _DELTA_Y * 0.5);
        return temperature;
    };

    inline float concentrationDiffInterpolated(float x, float y) const {
        float concentration = concentration_diff.valueInterpolated(x - _DELTA_X * 0.5, y - _DELTA_Y * 0.5);
        return concentration;
    };

    inline CellType cellType(int i, int j) const {
        i = clamp(i, 0, _SIZE_X - 1);
        j = clamp(j, 0, _SIZE_Y - 1);
        return _cell_type_buffer.value(i, j);
    };

    inline float divVelX(int i, int j) const {
        return (_vel_x_front_buffer->value(i + 1, j) - _vel_x_front_buffer->value(i, j)) / _DELTA_X;
    };

    inline float divVelY(int i, int j) const {
        return (_vel_y_front_buffer->value(i, j + 1) - _vel_y_front_buffer->value(i, j)) / _DELTA_Y;
    };

    // glm::dmat2 computeVelocityGradientMatrix(int i, int j);

    // Setters
    inline void setTemperature(int i, int j, float temperature) { (*_temperature_front_buffer)(i, j) = temperature; };

    inline void setTemperatureBackBuffer(int i, int j, float temperature) {
        (*_temperature_back_buffer)(i, j) = temperature;
    };

    inline void setConcentration(int i, int j, float concentration) {
        (*_concentration_front_buffer)(i, j) = concentration;
    };

    inline void setConcentrationBackBuffer(int i, int j, float concentration) {
        (*_concentration_back_buffer)(i, j) = concentration;
    };

    inline void setVelXHalfIndexed(int i, int j, float vel_x) { (*_vel_x_front_buffer)(i, j) = vel_x; };

    inline void setVelYHalfIndexed(int i, int j, float vel_y) { (*_vel_y_front_buffer)(i, j) = vel_y; };

    inline void setVelXBackBuffer(int i, int j, float vel_x) {
        (*_vel_x_back_buffer)(i, j) = vel_x;
        (*_vel_x_back_buffer)(i + 1, j) = vel_x;
    };

    inline void setVelYBackBuffer(int i, int j, float vel_y) {
        (*_vel_y_back_buffer)(i, j) = vel_y;
        (*_vel_y_back_buffer)(i, j + 1) = vel_y;
    };

    inline void setVelXBackBufferHalfIndexed(int i, int j, float vel_x) { (*_vel_x_back_buffer)(i, j) = vel_x; };

    inline void setVelYBackBufferHalfIndexed(int i, int j, float vel_y) { (*_vel_y_back_buffer)(i, j) = vel_y; };

    inline void setCellType(int i, int j, CellType cell_type) { _cell_type_buffer(i, j) = cell_type; };

    inline void addToVelXInterpolated(float x, float y, float vel_x) {
        _vel_x_back_buffer->addToValueInterpolated(x, y - 0.5 * _DELTA_Y, vel_x);
        // -0.5 Due to the MAC grid structure
    };

    inline void addToVelYInterpolated(float x, float y, float vel_y) {
        _vel_y_back_buffer->addToValueInterpolated(x - 0.5 * _DELTA_X, y, vel_y);
        // -0.5 Due to the MAC grid structure
    };

    inline void addToTemperatureInterpolated(float x, float y, float temperature) {
        _temperature_back_buffer->addToValueInterpolated(x, y, temperature);
    };

    inline void addToConcentrationInterpolated(float x, float y, float concentration) {
        _concentration_back_buffer->addToValueInterpolated(x, y, concentration);
    };

    inline void setVelXInterpolated(float x, float y, float vel_x) {
        _vel_x_back_buffer->setValueInterpolated(x, y - 0.5 * _DELTA_Y, vel_x);
        // -0.5 Due to the MAC grid structure
    };

    inline void setVelYInterpolated(float x, float y, float vel_y) {
        _vel_y_back_buffer->setValueInterpolated(x - 0.5 * _DELTA_X, y, vel_y);
        // -0.5 Due to the MAC grid structure
    };

    inline void setTemperatureInterpolated(float x, float y, float temperature) {
        _temperature_back_buffer->setValueInterpolated(x, y, temperature);
    };

    inline void setConcentrationInterpolated(float x, float y, float concentration) {
        _concentration_back_buffer->setValueInterpolated(x, y, concentration);
    };

    void swapVelocityBuffers();
    void swapTemperatureBuffers();
    void swapConcentrationBuffers();

    private:
    // Velocity grids need front buffers and back buffers
    // Normally read from front buffers and write to back buffers, then swap
    std::unique_ptr< Grid<float> > _vel_x_front_buffer;
    std::unique_ptr< Grid<float> > _vel_y_front_buffer;

    std::unique_ptr< Grid<float> > _vel_x_back_buffer;
    std::unique_ptr< Grid<float> > _vel_y_back_buffer;

    std::unique_ptr< Grid<float> > _temperature_front_buffer;
    std::unique_ptr< Grid<float> > _temperature_back_buffer;

    std::unique_ptr< Grid<float> > _concentration_front_buffer;
    std::unique_ptr< Grid<float> > _concentration_back_buffer;

    Grid<float> _vel_x_previous;
    Grid<float> _vel_y_previous;
    Grid<float> temperature_previous;
    Grid<float> concentration_previous;


    Grid<float> _vel_x_diff;
    Grid<float> _vel_y_diff;
    Grid<float> temperature_diff;
    Grid<float> concentration_diff;

    Grid<CellType> _cell_type_buffer;
};

#endif
