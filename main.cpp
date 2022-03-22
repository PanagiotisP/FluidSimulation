// 2D Fluid Simulation for now
// No rectangular bounding box, 4 free surfaces for now
#include "FluidSimulator.h"
#include "array2.h"
#include "array2_utils.h"
#include "vec.h"

#include <SFML/Graphics.hpp>
#ifdef _WIN32
    #include <Windows.h>
#else
    #include <unistd.h>
#endif

// #define EULER_GRID
#ifndef EULER_GRID
    #define PIC
#endif
#define SHOW_PARTICLES true
#define SHOW_VECTORS   false
#define INIT_TYPE      InitialisationType::VortexSmoke
#define PI             3.14159265f

// flip method
//� Transfer particle values qp to the grid qi, j, k, through equations like(7.1)
// or (7.2), and extrapolate on the grid as necessary.
//� Save the grid values qi, j, k.
//� Compute all other terms on the grid, such as pressure projection, to
// get an updated qnew
// i, j, k.
//� For each particle, interpolate the change qi, j, k = qnew
// i, j, k ? qi, j, k from
// the grid to add it to the particle�s value.
//� Advect the particles in the grid velocity field.

enum class InitialisationType { CenterSmoke, VortexSmoke, CollidingSmoke };

// Time span of each frame
float t_frame = 0.01f;

int max_frames_n = 100;

// Grid dimensions
int nx = 20;
int ny = 20;
float width = nx;

void euler_grid_initialisation(MacGrid &grid, std::vector<ScalarSource> &tempSrcs,
                               std::vector<ScalarSource> &concentrationSrcs, std::vector<VectorSource> &velocitySrcs,
                               InitialisationType type) {
    Vec2f center_point(((grid.sizeX() - 1) / 2 + 0.5) * grid.deltaX(), ((grid.sizeY() - 1) / 2 + 0.5) * grid.deltaY());

    Vec2f top_right(((grid.sizeX() - 2) + 0.5) * grid.deltaX(), ((grid.sizeY() - 2) + 0.5) * grid.deltaY());
    Vec2f bot_right(((grid.sizeX() - 2) + 0.5) * grid.deltaX(), (1 + 0.5) * grid.deltaY());
    Vec2f bot_left((1 + 0.5) * grid.deltaX(), (1 + 0.5) * grid.deltaY());
    Vec2f top_left((1 + 0.5) * grid.deltaX(), ((grid.sizeY() - 2) + 0.5) * grid.deltaY());

    Vec2f middle_top(((grid.sizeX() - 1) / 2 + 0.5) * grid.deltaX(), ((grid.sizeY() - 2) + 0.5) * grid.deltaY());
    Vec2f middle_right(((grid.sizeX() - 2) + 0.5) * grid.deltaX(), ((grid.sizeY() - 1) / 2 + 0.5) * grid.deltaY());
    Vec2f middle_bot(((grid.sizeX() - 1) / 2 + 0.5) * grid.deltaX(), (1 + 0.5) * grid.deltaY());
    Vec2f middle_left((1 + 0.5) * grid.deltaX(), ((grid.sizeY() - 1) / 2 + 0.5) * grid.deltaY());

    VectorSource vel_src1;
    VectorSource vel_src2;
    VectorSource vel_src3;
    VectorSource vel_src4;
    ScalarSource conc_src1;
    ScalarSource conc_src2;
    ScalarSource conc_src3;
    ScalarSource conc_src4;
    float velocity;
    float concentration;
    switch (type) {
        case InitialisationType::CenterSmoke:
            velocity = 200;
            concentration = 50;
            tempSrcs = std::vector<ScalarSource>(1, { center_point[0], center_point[1], ZERO_CELCIUS + 30 });
            concentrationSrcs = std::vector<ScalarSource>(1, { center_point[0], center_point[1], concentration });
            velocitySrcs = std::vector<VectorSource>(1, { center_point[0], center_point[1], velocity, 0 });
            return;
        case InitialisationType::VortexSmoke:
            velocity = 500;
            concentration = 70;

            vel_src1 = { top_right[0], top_right[1], 0, -velocity };
            vel_src2 = { bot_right[0], bot_right[1], -velocity, 0 };
            vel_src3 = { bot_left[0], bot_left[1], 0, velocity };
            vel_src4 = { top_left[0], top_left[1], velocity, 0 };

            conc_src1 = { top_right[0], top_right[1], concentration };
            conc_src2 = { bot_right[0], bot_right[1], concentration };
            conc_src3 = { bot_left[0], bot_left[1], concentration };
            conc_src4 = { top_left[0], top_left[1], concentration };

            concentrationSrcs = std::vector<ScalarSource> { conc_src1, conc_src2, conc_src3, conc_src4 };
            velocitySrcs = std::vector<VectorSource> { vel_src1, vel_src2, vel_src3, vel_src4 };
            return;
        case InitialisationType::CollidingSmoke:
            velocity = 100;
            concentration = 50;

            vel_src1 = { middle_top[0], middle_top[1], 0, -velocity };
            vel_src2 = { middle_bot[0], middle_bot[1], 0, velocity };
            vel_src3 = { middle_left[0], middle_left[1], velocity, 0 };
            vel_src4 = { middle_right[0], middle_right[1], -velocity, 0 };

            conc_src1 = { middle_top[0], middle_top[1], concentration };
            conc_src2 = { middle_bot[0], middle_bot[1], concentration };
            conc_src3 = { middle_left[0], middle_left[1], concentration };
            conc_src4 = { middle_right[0], middle_right[1], concentration };

            concentrationSrcs = std::vector<ScalarSource> { // conc_src1,
                                                            // conc_src2,
                                                            conc_src3, conc_src4
            };
            velocitySrcs = std::vector<VectorSource> { // vel_src1,
                                                       // vel_src2,
                                                       vel_src3, vel_src4
            };
            return;
    }
}

void pic_initialisation(FluidDomain &domain, InitialisationType type) {
    MacGrid &grid = domain.grid();
    Vec2f center_point(grid.lengthX() / 2, grid.lengthY() / 2);

    Vec2f top_right(grid.lengthX()-(1 + 0.5) * grid.deltaX(), grid.lengthY() - (1 + 0.5) * grid.deltaY());
    Vec2f bot_right(grid.lengthX()-(1 + 0.5) * grid.deltaX(), (1 + 0.5) * grid.deltaY());
    Vec2f bot_left((1 + 0.5) * grid.deltaX(), (1 + 0.5) * grid.deltaY());
    Vec2f top_left((1 + 0.5) * grid.deltaX(), grid.lengthY() - (1 + 0.5) * grid.deltaY());

    Vec2f middle_top(grid.lengthX() / 2, grid.lengthY() - (1 + 0.5) * grid.deltaY());
    Vec2f middle_right(grid.lengthX() -(1 + 0.5) *  grid.deltaX(), grid.lengthY() / 2);
    Vec2f middle_bot(grid.lengthX() / 2, (1 + 0.5) * grid.deltaY());
    Vec2f middle_left((1 + 0.5) * grid.deltaX(), grid.lengthY() / 2);

    float velocity;
    float concentration;
    // particles per grid cell per second
    int particle_generation_rate = 400;

    FluidSource src1;
    FluidSource src2;
    FluidSource src3;
    FluidSource src4;
    switch (type) {
        case InitialisationType::CenterSmoke:
            velocity = 200;
            concentration = 90;
            src1 = FluidSource(velocity, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                                  BBox<float> { center_point[0], center_point[0], center_point[1], center_point[1] },
                                  FluidDomain::ambient_temp, concentration);
            domain.addFluidSource(src1);
            return;
        case InitialisationType::VortexSmoke: 
            velocity = 500;
            concentration = 70;

            src1 = FluidSource(velocity, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                               BBox<float> { top_left[0], top_left[0], top_left[1], top_left[1] },
                               FluidDomain::ambient_temp, concentration);

            src2 = FluidSource(0.f, velocity, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                               BBox<float> { bot_left[0], bot_left[0], bot_left[1], bot_left[1] },
                               FluidDomain::ambient_temp, concentration);

            src3 = FluidSource(-velocity, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                               BBox<float> { bot_right[0], bot_right[0], bot_right[1], bot_right[1] },
                               FluidDomain::ambient_temp, concentration);

            src4 = FluidSource(0.f, -velocity, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                               BBox<float> { top_right[0], top_right[0], top_right[1], top_right[1] },
                               FluidDomain::ambient_temp, concentration);

            domain.addFluidSource(src1);
            domain.addFluidSource(src2);
            domain.addFluidSource(src3);
            domain.addFluidSource(src4);
            return;
        case InitialisationType::CollidingSmoke:
            velocity = 500;
            concentration = 125;

            src1 = FluidSource(velocity, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                               BBox<float> { middle_left[0], middle_left[0], middle_left[1], middle_left[1] },
                               FluidDomain::ambient_temp, concentration);

            src2 = FluidSource(0.f, velocity, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                               BBox<float> { middle_bot[0], middle_bot[0], middle_bot[1], middle_bot[1] },
                               FluidDomain::ambient_temp, concentration);

            src3 = FluidSource(-velocity, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                               BBox<float> { middle_right[0], middle_right[0], middle_right[1], middle_right[1] },
                               FluidDomain::ambient_temp, concentration);

            src4 = FluidSource(0.f, -velocity, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                               BBox<float> { middle_top[0], middle_top[0], middle_top[1], middle_top[1] },
                               FluidDomain::ambient_temp, concentration);

            domain.addFluidSource(src1);
            domain.addFluidSource(src2);
            domain.addFluidSource(src3);
            domain.addFluidSource(src4);
            return;
    }
}

int main(int argc, char **argv) {
    int frame = 0;
    float grid_square_width = nx * 10;

    FluidSimulator sim;
    FluidDomain fluidDomain(nx, ny, width, width);

    std::vector<ScalarSource> tempSrcs;
    std::vector<ScalarSource> concentrationSrcs;
    std::vector<VectorSource> velocitySrcs;
    bool isTempSrcsActive = false;
    bool isConcentrationSrcsActive = false;
    bool isVelocitySrcsActive = false;

    InitialisationType initType = INIT_TYPE;
#if defined EULER_GRID
    euler_grid_initialisation(fluidDomain.grid(), tempSrcs, concentrationSrcs, velocitySrcs, initType);
#elif defined PIC
    pic_initialisation(fluidDomain, initType);
#endif

    sf::RenderWindow window(sf::VideoMode(800, 800), "SFML works!");

    int grid_width_x, grid_width_y;
    fluidDomain.grid().sizeX() > fluidDomain.grid().sizeY() ? grid_width_x = 540 : grid_width_y = 540;
    fluidDomain.grid().sizeX() > fluidDomain.grid().sizeY()
        ? grid_width_y = 540 * fluidDomain.grid().sizeY() / fluidDomain.grid().sizeX()
        : grid_width_x = 540 * fluidDomain.grid().sizeX() / fluidDomain.grid().sizeY();
    assert(grid_width_x < 800);


    Array2<sf::RectangleShape> render_grid(fluidDomain.grid().sizeX(), fluidDomain.grid().sizeY());
    Array2<sf::RectangleShape> vectors(fluidDomain.grid().sizeX(), fluidDomain.grid().sizeY());
    std::vector<sf::RectangleShape> render_particles;

    sf::RectangleShape grid_rectangle(
        sf::Vector2f(grid_width_x / fluidDomain.grid().sizeX(), grid_width_y / fluidDomain.grid().sizeY()));
    grid_rectangle.setOutlineThickness(0);
    grid_rectangle.setOutlineColor(sf::Color::Green);

    render_grid.assign(grid_rectangle);


    sf::RectangleShape vector_arrow(sf::Vector2f(0, 2));
    vector_arrow.setFillColor(sf::Color::White);

    vectors.assign(vector_arrow);

    sf::Vector2f offset((800 - grid_width_x) / 2, (800 - grid_width_y) / 2);
    for (int j = 0; j < fluidDomain.grid().sizeY(); ++j)
        for (int i = 0; i < fluidDomain.grid().sizeX(); ++i) {
            if (fluidDomain.grid().cellType(i, j) == SOLID)
                render_grid(i, j).setFillColor(sf::Color::Yellow);
            else if (fluidDomain.grid().cellType(i, j) == LIQUID)
                render_grid(i, j).setFillColor(sf::Color::Black);
            else
                render_grid(i, j).setFillColor(sf::Color::Magenta);
            render_grid(i, j).setPosition(offset.x + i * grid_width_x / fluidDomain.grid().sizeX(),
                                          offset.y + grid_width_y
                                              - (j + 1) * grid_width_y / fluidDomain.grid().sizeY());
            vectors(i, j).setPosition(offset.x + (i + 0.5) * grid_width_x / fluidDomain.grid().sizeX(),
                                      offset.y + grid_width_y - (j + 0.5) * grid_width_y / fluidDomain.grid().sizeY());
        }


    while (window.isOpen() && (true || frame < 50)) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) window.close();
        }

        for (int j = 0; j < fluidDomain.grid().sizeY(); ++j) {
            for (int i = 0; i < fluidDomain.grid().sizeX(); ++i) {
                Vec2f point((i + 0.5) * fluidDomain.grid().deltaX(), (j + 0.5) * fluidDomain.grid().deltaY());
                Vec2f velocity_vector =
                    Vec2f(fluidDomain.grid().velXInterpolated(point[0], point[1] + 0.5 * fluidDomain.grid().deltaY()),
                          fluidDomain.grid().velYInterpolated(point[0] + 0.5 * fluidDomain.grid().deltaX(), point[1]));

                if (fluidDomain.grid().cellType(i, j) == LIQUID)
                    render_grid(i, j).setFillColor(sf::Color(255, 255, 255,
                                                             255
                                                                 * fluidDomain.grid().concentrationInterpolated(
                                                                     point[0] ,/*+ 0.5 * fluidDomain.grid().deltaX()*/
                                                                     point[1] /*+ 0.5 * fluidDomain.grid().deltaY()*/)));

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
                if (SHOW_VECTORS) {
                    vectors(i, j).setRotation(angle);
                    vectors(i, j).setSize(sf::Vector2f(mag(velocity_vector), 2));
                }
            }
        }
        render_particles.clear();
        for (auto it = fluidDomain.particleSet().begin(); it != fluidDomain.particleSet().end(); ++it) {
            sf::RectangleShape particle(sf::Vector2f(3, 3));
            particle.setFillColor(sf::Color::Magenta);
            particle.setPosition(
                offset.x + (*it)->posX() * grid_width_x / fluidDomain.grid().sizeX() / fluidDomain.grid().deltaX(),
                offset.y + grid_width_y
                    - (*it)->posY() * grid_width_y / fluidDomain.grid().sizeY() / fluidDomain.grid().deltaY());
            render_particles.push_back(particle);
        }
        window.clear(sf::Color::Black);
        for (int i = 0; i < render_grid.a.size(); ++i) {
            assert(render_grid.a.size() == vectors.a.size());
            window.draw(render_grid.a[i]);
            if (SHOW_VECTORS) window.draw(vectors.a[i]);
        }

        for (auto it = render_particles.begin(); it != render_particles.end(); ++it) {
            if (SHOW_PARTICLES) window.draw(*it);
        }

        window.display();

        if (frame % 200 == 0) {
            isVelocitySrcsActive = !isVelocitySrcsActive;
            isConcentrationSrcsActive = !isConcentrationSrcsActive;
        }
#if defined EULER_GRID
        sim.advance_eulerian_grid(fluidDomain, t_frame, &tempSrcs, &concentrationSrcs, &velocitySrcs, isTempSrcsActive,
                                  isConcentrationSrcsActive, isVelocitySrcsActive);
#elif defined PIC
        sim.advance_flip_pic(fluidDomain, t_frame, 0.98);
#endif
        // Sleep(50);
        frame++;
    }

    return 0;
}