// 2D Fluid Simulation for now
// No rectangular bounding box, 4 free surfaces for now
#include "FluidDomain.h"
#include "FluidSimulator.h"
#include "Renderer.h"
#include "array2.h"
#include "array2_utils.h"
#include "render_utils.h"
#include "vec.h"

#include <SFML/Graphics.hpp>
#include <array>
#include <iomanip>
#ifdef _WIN32
    #include <Windows.h>
#else
    #include <unistd.h>
#endif

// #define EULER_GRID
#ifndef EULER_GRID
    #define PIC_FLIP
#endif

#define INIT_TYPE InitialisationType::CenterSmoke

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

// Time span of each frame
float t_frame = 0.016f;

int max_frames_n = 100;

// Grid dimensions
int nx = 64;
int ny = 64;
float widthY = 1;
float widthX = 1;

enum class InitialisationType { CenterSmoke, VortexSmoke, CollidingSmoke, WaterFall };

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

    Vec2f top_right(grid.lengthX() - (1 + 0.5) * grid.deltaX(), grid.lengthY() - (1 + 0.5) * grid.deltaY());
    Vec2f bot_right(grid.lengthX() - (1 + 0.5) * grid.deltaX(), (1 + 0.5) * grid.deltaY());
    Vec2f bot_left(0.5 * grid.deltaX(), (1 + 0.5) * grid.deltaY());
    Vec2f top_left(0.5 * grid.deltaX(), grid.lengthY() - (1 + 0.5) * grid.deltaY());

    Vec2f middle_top((grid.sizeX() / 2 + 0.5) * grid.deltaX(), grid.lengthY() - 2 * grid.deltaY());
    Vec2f middle_right(grid.lengthX() - 2.5 * grid.deltaX(), (grid.sizeY() / 2 + 0.5) * grid.deltaY());
    Vec2f middle_bot((grid.sizeX() / 2 + 0.5) * grid.deltaX(), 2 * grid.deltaY());
    Vec2f middle_left(2.5 * grid.deltaX(), (grid.sizeY() / 2 + 0.5) * grid.deltaY());

    float velocity;
    float concentration;
    // particles per grid cell per second
    int particle_generation_rate = 625;

    FluidSource src1(grid.sizeX(), grid.sizeY(), grid.lengthX(), grid.lengthY());
    FluidSource src2(grid.sizeX(), grid.sizeY(), grid.lengthX(), grid.lengthY());
    FluidSource src3(grid.sizeX(), grid.sizeY(), grid.lengthX(), grid.lengthY());
    FluidSource src4(grid.sizeX(), grid.sizeY(), grid.lengthX(), grid.lengthY());
    LevelSet l(grid.sizeX(), grid.sizeY(), grid.lengthX(), grid.lengthY());
    switch (type) {
        case InitialisationType::CenterSmoke:
            velocity = -10;
            concentration = 50;

            l.construct_from_known_geometry(distance_from_sphere(center_point[0], center_point[1], 0.1f));
            //  distance_from_axis_aligned_box(
            //     BBox<float> { center_point[0] - 6.0f * grid.deltaX(), center_point[0] + 6.0f * grid.deltaX(),
            //                   center_point[1] - 6.0f * grid.deltaY(), center_point[1] + 6.0f * grid.deltaY() }));
            src1 = FluidSource(l, velocity, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                               FluidDomain::ambient_temp, concentration);
            domain.addFluidSource(src1);
            return;
        case InitialisationType::VortexSmoke:
            velocity = 100;
            concentration = 50;

            // src1 = FluidSource(velocity, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
            //                    BBox<float> { top_left[0], top_left[0], top_left[1], top_left[1] },
            //                    FluidDomain::ambient_temp, concentration);

            // src2 = FluidSource(0.f, velocity, grid.deltaX(), grid.deltaY(), particle_generation_rate,
            //                    BBox<float> { bot_left[0], bot_left[0], bot_left[1], bot_left[1] },
            //                    FluidDomain::ambient_temp, concentration);

            // src3 = FluidSource(-velocity, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
            //                    BBox<float> { bot_right[0], bot_right[0], bot_right[1], bot_right[1] },
            //                    FluidDomain::ambient_temp, concentration);

            // src4 = FluidSource(0.f, -velocity, grid.deltaX(), grid.deltaY(), particle_generation_rate,
            //                    BBox<float> { top_right[0], top_right[0], top_right[1], top_right[1] },
            //                    FluidDomain::ambient_temp, concentration);

            // domain.addFluidSource(src1);
            // domain.addFluidSource(src2);
            // domain.addFluidSource(src3);
            // domain.addFluidSource(src4);
            return;
        case InitialisationType::CollidingSmoke:
            velocity = 500;
            concentration = 125;

            // src1 = FluidSource(velocity, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
            //                    BBox<float> { middle_left[0], middle_left[0], middle_left[1], middle_left[1] },
            //                    FluidDomain::ambient_temp, concentration);

            // src2 = FluidSource(0.f, velocity, grid.deltaX(), grid.deltaY(), particle_generation_rate,
            //                    BBox<float> { middle_bot[0], middle_bot[0], middle_bot[1], middle_bot[1] },
            //                    FluidDomain::ambient_temp, concentration);

            // src3 = FluidSource(-velocity, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
            //                    BBox<float> { middle_right[0], middle_right[0], middle_right[1], middle_right[1] },
            //                    FluidDomain::ambient_temp, concentration);

            // src4 = FluidSource(0.f, 0, grid.deltaX(), grid.deltaY(), particle_generation_rate,
            //                    BBox<float> { middle_top[0], middle_top[0], middle_top[1], middle_top[1] },
            //                    FluidDomain::ambient_temp, concentration);

            // domain.addFluidSource(src1);
            // domain.addFluidSource(src2);
            // domain.addFluidSource(src3);
            // domain.addFluidSource(src4);
            return;
        case InitialisationType::WaterFall:
            velocity = 0;
            concentration = 70;

            l.construct_from_known_geometry(distance_from_axis_aligned_box(
                BBox<float> { middle_bot[0] - 30.f * grid.deltaX(), middle_bot[0] + 30.f * grid.deltaX(),
                              middle_bot[1] - 1.0f * grid.deltaY(), middle_bot[1] + 1.0f * grid.deltaY() }));
            src1 = FluidSource(l, velocity, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                               FluidDomain::ambient_temp, concentration);

            LevelSet l2(l);
            l2.construct_from_known_geometry(distance_from_axis_aligned_box(
                BBox<float> { middle_top[0] - 5.0f * grid.deltaX(), middle_top[0] + 5.0f * grid.deltaX(),
                              middle_top[1] - 5.0f * grid.deltaY(), middle_top[1] + 5.0f * grid.deltaY() }));
            src2 = FluidSource(l2, 0.f, -10.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
                               FluidDomain::ambient_temp, concentration);

            // src3 = FluidSource(0, 0.f, grid.deltaX(), grid.deltaY(), particle_generation_rate,
            //                    BBox<float> { top_right[0], top_right[0], top_right[1], top_right[1] },
            //                    FluidDomain::ambient_temp, concentration);

            // src4 = FluidSource(0, 0, grid.deltaX(), grid.deltaY(), particle_generation_rate,
            //                    BBox<float> { top_left[0], top_left[0], top_left[1], top_left[1] },
            //                    FluidDomain::ambient_temp, concentration);

            domain.addFluidSource(src1);
            domain.addFluidSource(src2);
            // domain.addFluidSource(src3);
            // domain.addFluidSource(src4);
            return;
    }
}

int main(int argc, char **argv) {
    std::cout.precision(2);

    int frame = 0;

    LevelSet solid_level_set(nx, ny, widthX, widthY);
    solid_level_set.construct_from_known_geometry(distance_from_axis_aligned_box(
        BBox<float>(0.5f * solid_level_set.deltaX(),
                    (solid_level_set.sizeX() - 2) * solid_level_set.deltaX() + 0.5f * solid_level_set.deltaX(),
                    0.5f * solid_level_set.deltaY(),
                    (solid_level_set.sizeY() - 2) * solid_level_set.deltaY() + 0.5f * solid_level_set.deltaY()),
        true));
    // solid_level_set.construct_from_edges(std::vector({ Vec2f(0.5, 0.5), Vec2f(9.5, 0.5), Vec2f(9.5, 9.5) }), true);
    // solid_level_set.construct_from_known_geometry(distance_from_sphere(widthX/2, widthY/2, widthX/2, true));
    FluidSimulator sim;
    FluidDomain fluidDomain(nx, ny, widthX, widthY);
    fluidDomain.setSolidLevelSet(solid_level_set);

    Renderer renderer(nx, ny);
    renderer.initialisation();

    std::vector<ScalarSource> tempSrcs;
    std::vector<ScalarSource> concentrationSrcs;
    std::vector<VectorSource> velocitySrcs;
    bool isTempSrcsActive = false;
    bool isConcentrationSrcsActive = false;
    bool isVelocitySrcsActive = false;

    InitialisationType initType = INIT_TYPE;
#if defined EULER_GRID
    euler_grid_initialisation(fluidDomain.grid(), tempSrcs, concentrationSrcs, velocitySrcs, initType);
#elif defined PIC_FLIP
    pic_initialisation(fluidDomain, initType);
#endif


    sf::RenderWindow window = sf::RenderWindow(sf::VideoMode(WINDOW_SIZE, WINDOW_SIZE), "SFML works!");
    while (window.isOpen() && (true || frame < 50)) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) window.close();
        }

        window.clear(sf::Color::Black);
        renderer.draw(fluidDomain, window);
        window.display();

        if (frame % 2000 == 100) {
            // isVelocitySrcsActive = !isVelocitySrcsActive;
            // isConcentrationSrcsActive = !isConcentrationSrcsActive;
            fluidDomain.removeFluidSource(0);
        }
        if (frame % 2000 == 1) { fluidDomain.removeFluidSource(0); }

#if defined EULER_GRID
        sim.advance_eulerian_grid(fluidDomain, t_frame, &tempSrcs, &concentrationSrcs, &velocitySrcs, isTempSrcsActive,
                                  isConcentrationSrcsActive, isVelocitySrcsActive);
#elif defined PIC_FLIP
        sim.advance_flip_pic(fluidDomain, t_frame, 0.95);
#endif
        // Sleep(50);
        frame++;
    }

    return 0;
}