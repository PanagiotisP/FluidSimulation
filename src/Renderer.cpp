#include "Renderer.h"

#include <openvdb/tools/Interpolation.h>

Renderer::Renderer(int sim_size_x, int sim_size_y):
 sim_size_x(sim_size_x), sim_size_y(sim_size_y), render_grid(sim_size_x, sim_size_x), vectors(sim_size_y, sim_size_y) {
    sim_size_x > sim_size_y ? grid_width_x = 540 : grid_width_y = 540;
    sim_size_x > sim_size_y ? grid_width_y = 540 *sim_size_y / sim_size_x
                            : grid_width_x = 540 * sim_size_x / sim_size_y;
    offset = sf::Vector2f((WINDOW_SIZE - grid_width_x) / 2, (WINDOW_SIZE - grid_width_y) / 2);
}

void Renderer::initialisation() {
    float grid_cell_width_x = grid_width_x / sim_size_x;
    float grid_cell_width_y = grid_width_y / sim_size_y;
    sf::RectangleShape grid_rectangle(sf::Vector2f(grid_cell_width_x, grid_cell_width_y));
    grid_rectangle.setOutlineColor(sf::Color::Green);
    grid_rectangle.setOutlineThickness(0);
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
    MacGrid::BoxSampler vel_sampler(fluidDomain.grid().velFront()->getAccessor(),
                                    fluidDomain.grid().velFront()->transform());
    if (SHOW_VECTORS) {
        for (int j = 0; j < sim_size_y; ++j) {
            for (int i = 0; i < sim_size_x; ++i) {
                Vec2f point((i + 0.5) * fluidDomain.voxelSize(), (j + 0.5) * fluidDomain.voxelSize());
                Vec2f velocity_vector =
                    Vec2f(fluidDomain.grid().velInterpolatedW(vel_sampler, sim_size_x / 2, point[1], point[0])[2],
                          fluidDomain.grid().velInterpolatedW(vel_sampler, sim_size_x / 2, point[1], point[0])[1]);

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
        auto fluidAccessor = fluidDomain.fluidLevelSet().getAccessor();
        auto solidAccessor = fluidDomain.solidLevelSet().getAccessor();
        LevelSet::BoxSampler fluidSampler(fluidAccessor, fluidDomain.fluidLevelSet().getLevelSet()->transform());
        LevelSet::BoxSampler solidSampler(solidAccessor, fluidDomain.solidLevelSet().getLevelSet()->transform());

        for (int j = 0; j < sim_size_y; ++j) {
            for (int i = 0; i < sim_size_x; ++i) {
                if (fluidAccessor.getValue(openvdb::Coord(sim_size_x / 2, j, i)) < 0) {
#if defined SMOKE
                    render_grid(i, j).setFillColor(
                        sf::Color(255, 255, 255,
                                  255
                                      * fluidAccessor.getValue(openvdb::Coord(sim_size_x / 2, j, i))
                                            concentrationInterpolated(point[0], point[1])));
#else
                    render_grid(i, j).setFillColor(sf::Color::Cyan);
#endif
                } else if (fluidAccessor.getValue(openvdb::Coord(sim_size_x / 2, j, i)) > 0
                           && solidAccessor.getValue(openvdb::Coord(sim_size_x / 2, j, i)) > 0) { // Check for AIR cell
                    render_grid(i, j).setFillColor(sf::Color::White);
                } else {
                    render_grid(i, j).setFillColor(sf::Color::Yellow);
                }
                window.draw(render_grid(i, j));
            }
        }
    }

    if (SHOW_PARTICLES) {
        for (auto it = fluidDomain.particleSet().begin(); it != fluidDomain.particleSet().end(); ++it) {
            // if (2 < (*it).p[0] && (*it).p[0] < 3) {
            sf::CircleShape particle(2); // grid_width_x / 2 / sim_size_x);
            particle.setFillColor(sf::Color::Blue);
            particle.setPosition(offset.x + it->pos()[2] * grid_width_x / sim_size_x,
                                 offset.y + grid_width_y - it->pos()[1] * grid_width_y / sim_size_y);
            window.draw(particle);
            // }
        }
    }

    if (SHOW_LEVEL_SETS) {
        // auto solid_vertices = marching_cubes_2d(fluidDomain.solidLevelSet(), grid_width_x /
        // sim_size_x,
        //                                         grid_width_y / fluidDomain.grid().sizeY(), offset.x,
        //                                         offset.y + grid_width_y, sf::Color::Yellow);

        // auto fluid_vertices = marching_cubes_2d(fluidDomain.fluidLevelSet(), grid_width_x /
        // sim_size_x,
        //                                         grid_width_y / fluidDomain.grid().sizeY(), offset.x,
        //                                         offset.y + grid_width_y, sf::Color::Blue);

        // for (auto it = solid_vertices.begin(); it != solid_vertices.end(); ++it) {
        //     sf::Vertex line[2] = { (*it)[0], (*it)[1] };
        //     window.draw(line, 2, sf::Lines);
        // }

        // for (auto it = fluid_vertices.begin(); it != fluid_vertices.end(); ++it) {
        //     sf::Vertex line[2] = { (*it)[0], (*it)[1] };
        //     window.draw(line, 2, sf::Lines);
        // }
    }
}