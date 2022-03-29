#include "Renderer.h"

Renderer::Renderer(int sim_size_x, int sim_size_y):
 sim_size_x(sim_size_x), sim_size_y(sim_size_y), render_grid(sim_size_x, sim_size_x), vectors(sim_size_y, sim_size_y){
    sim_size_x > sim_size_y ? grid_width_x = 540 : grid_width_y = 540;
    sim_size_x > sim_size_y ? grid_width_y = 540 *sim_size_y / sim_size_x
                            : grid_width_x = 540 * sim_size_x / sim_size_y;
    offset = sf::Vector2f((WINDOW_SIZE - grid_width_x) / 2, (WINDOW_SIZE - grid_width_y) / 2);
}

void Renderer::initialisation() {
    float grid_cell_width_x = grid_width_x / sim_size_x;
    float grid_cell_width_y = grid_width_y / sim_size_y;
    sf::RectangleShape grid_rectangle(sf::Vector2f(grid_cell_width_x, grid_cell_width_y));
    render_grid.assign(grid_rectangle);

    sf::RectangleShape vector_arrow(sf::Vector2f(0, 2));
    vector_arrow.setFillColor(sf::Color::Black);
    vectors.assign(vector_arrow);
    for (int j = 0; j < sim_size_y; ++j)
        for (int i = 0; i < sim_size_x; ++i) {
            render_grid(i, j).setPosition(offset.x + i * grid_cell_width_x,
                                          offset.y + grid_width_y - (j + 1) * grid_cell_width_y);

            // Inverse offset, so that the 1st row is rendered on the bottom
            vectors(i, j).setPosition(offset.x + (i + 0.5) * grid_cell_width_x,
                                      offset.y + grid_width_y - (j + 0.5) * grid_cell_width_y);
        }
}

void Renderer::draw(FluidDomain &fluidDomain, sf::RenderWindow &window) {
    if (SHOW_VECTORS) {
        for (int j = 0; j < sim_size_y; ++j) {
            for (int i = 0; i < sim_size_x; ++i) {
                Vec2f point((i + 0.5) * fluidDomain.grid().deltaX(), (j + 0.5) * fluidDomain.grid().deltaY());
                Vec2f velocity_vector = Vec2f(fluidDomain.grid().velXInterpolated(point[0], point[1]),
                                              fluidDomain.grid().velYInterpolated(point[0], point[1]));


                float angle = atan(fabs(velocity_vector[1]) / fabs(velocity_vector[0])) * 180.f / PI;

                if (velocity_vector[0] >= 0 && velocity_vector[1] >= 0) {
                    angle = -angle;
                } else if (velocity_vector[0] < 0 && velocity_vector[1] >= 0) {
                    angle = -(180 - angle);
                } else if (velocity_vector[0] < 0 && velocity_vector[1] < 0) {
                    angle = -(180 + angle);
                } else if (velocity_vector[0] >= 0 && velocity_vector[1] < 0) {
                    angle = angle;
                }
                vectors(i, j).setRotation(angle);
                vectors(i, j).setSize(sf::Vector2f(mag(velocity_vector), 2));

                window.draw(vectors(i, j));
            }
        }
    }

    if (SHOW_GRID_VOXELS) {
        for (int j = 0; j < sim_size_y; ++j) {
            for (int i = 0; i < sim_size_x; ++i) {
                if (fluidDomain.fluidLevelSet()(i, j) < 0) {
#if defined SMOKE
                    render_grid(i, j).setFillColor(sf::Color(
                        255, 255, 255, 255 * fluidDomain.grid().concentrationInterpolated(point[0], point[1])));
#else
                    render_grid(i, j).setFillColor(sf::Color::Cyan);
#endif
                } else if (fluidDomain.fluidLevelSet()(i, j) > 0
                           && fluidDomain.solidLevelSet()(i, j) > 0) { // Check for AIR cell
                    render_grid(i, j).setFillColor(sf::Color::White);
                }
                window.draw(render_grid(i, j));
            }
        }
    }

    if (SHOW_PARTICLES) {
        for (auto it = fluidDomain.particleSet().begin(); it != fluidDomain.particleSet().end(); ++it) {
            sf::CircleShape particle(2); // grid_width_x / 2 / sim_size_x);
            particle.setFillColor(sf::Color::Blue);
            particle.setPosition(offset.x + (*it)->posX() * grid_width_x / sim_size_x / fluidDomain.grid().deltaX(),
                                 offset.y + grid_width_y
                                     - (*it)->posY() * grid_width_y / sim_size_y / fluidDomain.grid().deltaY());
            window.draw(particle);
        }
    }

    if (SHOW_LEVEL_SETS) {
        auto solid_vertices = marching_cubes_2d(fluidDomain.solidLevelSet(), grid_width_x / fluidDomain.grid().sizeX(),
                                                grid_width_y / fluidDomain.grid().sizeY(), offset.x,
                                                offset.y + grid_width_y, sf::Color::Yellow);

        auto fluid_vertices = marching_cubes_2d(fluidDomain.fluidLevelSet(), grid_width_x / fluidDomain.grid().sizeX(),
                                                grid_width_y / fluidDomain.grid().sizeY(), offset.x,
                                                offset.y + grid_width_y, sf::Color::Blue);

        for (auto it = solid_vertices.begin(); it != solid_vertices.end(); ++it) {
            sf::Vertex line[2] = { (*it)[0], (*it)[1] };
            window.draw(line, 2, sf::Lines);
        }

        for (auto it = fluid_vertices.begin(); it != fluid_vertices.end(); ++it) {
            sf::Vertex line[2] = { (*it)[0], (*it)[1] };
            window.draw(line, 2, sf::Lines);
        }
    }
}

std::vector<std::array<sf::Vertex, 2>> Renderer::marching_cubes_2d(LevelSet &l, int grid_cell_width_y,
                                                                   int grid_cell_width_x, int offsetX, int offsetY,
                                                                   sf::Color colour) {
    std::vector<std::array<sf::Vertex, 2>> lines;
    for (int j = 0; j < l.sizeY() - 1; ++j) {
        for (int i = 0; i < l.sizeX() - 1; ++i) {
            float v00 = l(i, j);
            float v01 = l(i, j + 1);
            float v10 = l(i + 1, j);
            float v11 = l(i + 1, j + 1);

            // Build the binary number
            int case_num = (v00 > 0) * 1 + (v01 > 0) * 2 + (v10 > 0) * 4 + (v11 > 0) * 8;

            float frac_1;
            float frac_2;
            std::array<sf::Vertex, 2> line;
            switch (case_num) {
                case 0:
                case 15: break;
                // Single Corner
                case 1:
                case 14:
                    frac_1 = -v00 / (v10 - v00);
                    frac_2 = -v00 / (v01 - v00);
                    line[0] = sf::Vertex(sf::Vector2f(offsetX + (i + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + frac_2 + 0.5) * grid_cell_width_y),
                                         colour);
                    line[1] = sf::Vertex(sf::Vector2f(offsetX + (i + frac_1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + 0.5) * grid_cell_width_y),
                                         colour);
                    lines.push_back(line);
                    break;
                case 2:
                case 13:
                    frac_1 = -v01 / (v11 - v01);
                    frac_2 = -v00 / (v01 - v00);
                    line[0] = sf::Vertex(sf::Vector2f(offsetX + (i + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + frac_2 + 0.5) * grid_cell_width_y),
                                         colour);
                    line[1] = sf::Vertex(sf::Vector2f(offsetX + (i + frac_1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + 1 + 0.5) * grid_cell_width_y),
                                         colour);
                    lines.push_back(line);
                    break;
                case 4:
                case 11:
                    frac_1 = -v00 / (v10 - v00);
                    frac_2 = -v10 / (v11 - v10);
                    line[0] = sf::Vertex(sf::Vector2f(offsetX + (i + frac_1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + 0.5) * grid_cell_width_y),
                                         colour);
                    line[1] = sf::Vertex(sf::Vector2f(offsetX + (i + 1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + frac_2 + 0.5) * grid_cell_width_y),
                                         colour);
                    lines.push_back(line);
                    break;
                case 7:
                case 8:
                    frac_1 = -v01 / (v11 - v01);
                    frac_2 = -v10 / (v11 - v10);
                    line[0] = sf::Vertex(sf::Vector2f(offsetX + (i + 1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + frac_2 + 0.5) * grid_cell_width_y),
                                         colour);
                    line[1] = sf::Vertex(sf::Vector2f(offsetX + (i + frac_1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + 1 + 0.5) * grid_cell_width_y),
                                         colour);
                    lines.push_back(line);
                    break;
                    // Vertical line
                case 3:
                case 12:
                    frac_1 = -v00 / (v10 - v00);
                    frac_2 = -v01 / (v11 - v01);
                    line[0] = sf::Vertex(sf::Vector2f(offsetX + (i + frac_1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + 0.5) * grid_cell_width_y),
                                         colour);
                    line[1] = sf::Vertex(sf::Vector2f(offsetX + (i + frac_2 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + 1 + 0.5) * grid_cell_width_y),
                                         colour);
                    lines.push_back(line);
                    break;
                // Horizontal line
                case 5:
                case 10:
                    frac_1 = -v00 / (v01 - v00);
                    frac_2 = -v10 / (v11 - v10);
                    line[0] = sf::Vertex(sf::Vector2f(offsetX + (i + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + frac_1 + 0.5) * grid_cell_width_y),
                                         colour);
                    line[1] = sf::Vertex(sf::Vector2f(offsetX + (i + 1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + frac_2 + 0.5) * grid_cell_width_y),
                                         colour);
                    lines.push_back(line);
                    break;
                case 6:
                    frac_1 = -v01 / (v11 - v01);
                    frac_2 = -v00 / (v01 - v00);
                    line[0] = sf::Vertex(sf::Vector2f(offsetX + (i + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + frac_2 + 0.5) * grid_cell_width_y),
                                         colour);
                    line[1] = sf::Vertex(sf::Vector2f(offsetX + (i + frac_1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + 1 + 0.5) * grid_cell_width_y),
                                         colour);
                    lines.push_back(line);

                    frac_1 = -v00 / (v10 - v00);
                    frac_2 = -v10 / (v11 - v10);
                    line[0] = sf::Vertex(sf::Vector2f(offsetX + (i + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + frac_2 + 0.5) * grid_cell_width_y),
                                         colour);
                    line[1] = sf::Vertex(sf::Vector2f(offsetX + (i + frac_1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + 0.5) * grid_cell_width_y),
                                         colour);
                    lines.push_back(line);
                    break;
                case 9:
                    frac_1 = -v00 / (v10 - v00);
                    frac_2 = -v00 / (v01 - v00);
                    line[0] = sf::Vertex(sf::Vector2f(offsetX + (i + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + frac_2 + 0.5) * grid_cell_width_y),
                                         colour);
                    line[1] = sf::Vertex(sf::Vector2f(offsetX + (i + frac_1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + 0.5) * grid_cell_width_y),
                                         colour);
                    lines.push_back(line);
                    frac_1 = -v01 / (v11 - v01);
                    frac_2 = -v10 / (v11 - v10);
                    line[0] = sf::Vertex(sf::Vector2f(offsetX + (i + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + frac_2 + 0.5) * grid_cell_width_y),
                                         colour);
                    line[1] = sf::Vertex(sf::Vector2f(offsetX + (i + frac_1 + 0.5) * grid_cell_width_x,
                                                      offsetY - (j + 0.5) * grid_cell_width_y),
                                         colour);
                    lines.push_back(line);
                    break;
            }
        }
    }
    return lines;
}
