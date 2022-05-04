#pragma once

#include "FluidDomain.h"
#include "array2.h"
#include "render_utils.h"

#include <SFML/Graphics.hpp>

#define WINDOW_SIZE      800
#define SHOW_PARTICLES   true
#define SHOW_VECTORS     true
#define SHOW_GRID_VOXELS true
#define SHOW_LEVEL_SETS  false
// #define SMOKE
#define PI 3.14159265f

class Renderer {
public:
    Renderer(int sim_size_x, int sim_size_y);

    void initialisation();

    inline int widthX() { return grid_width_x; }
    inline int widthY() { return grid_width_y; }

    void draw(FluidDomain &fluidDomain, sf::RenderWindow &window);

private:

    std::vector<std::array<sf::Vertex, 2>> marching_cubes_2d(LevelSet &l, int grid_cell_width_y, int grid_cell_width_x,
                                                             int offsetX, int offsetY, sf::Color colour);

    int sim_size_x, sim_size_y;
    int grid_width_x, grid_width_y;
    Array2<sf::RectangleShape> render_grid;
    Array2<sf::RectangleShape> vectors;
    sf::Vector2f offset;
};