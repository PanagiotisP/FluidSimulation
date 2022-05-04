// 2D Fluid Simulation for now
// No rectangular bounding box, 4 free surfaces for now
#include "FluidDomain.h"
#include "FluidSimulator.h"
#include "Renderer.h"
#include "render_utils.h"
#include "vec.h"

#include <SFML/Graphics.hpp>
#include <array>
#include <filesystem>
#include <iomanip>
#include <openvdb/tools/FastSweeping.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/LevelSetPlatonic.h>
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/LevelSetUtil.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/tools/LevelSetMeasure.h>

#ifdef _WIN32
    #include <Windows.h>
#else
    #include <unistd.h>
#endif

// #define EULER_GRID
#ifndef EULER_GRID
    #define PIC_FLIP
#endif

#define INIT_TYPE InitialisationType::WaterFall


// using namespace openvdb;
// using namespace openvdb::tools;

// Time span of each frame
float t_frame = 0.008f;

int max_frames_n = 200;

// Grid dimensions
int nx = 64;
int ny = 64;
int nz = 64;
float width = 1;

enum class InitialisationType { CenterSmoke, VortexSmoke, CollidingSmoke, WaterFall };

void pic_initialisation(FluidDomain &domain, InitialisationType type, openvdb::math::Transform::Ptr i2w_transform) {
    MacGrid &grid = domain.grid();
    openvdb::Vec3f center_point(grid.sizeY() / 2, grid.sizeX() / 2, grid.sizeY() / 2);

    float velocity;
    float concentration;
    // particles per grid cell per second
    int particle_generation_rate = 625;

    FluidSource src1;
    FluidSource src2;
    FluidSource src3;
    FluidSource src4;
    openvdb::FloatGrid::Ptr l;
    switch (type) {
        case InitialisationType::CenterSmoke:
            velocity = 0;
            concentration = 50;

            src1 = FluidSource(l, openvdb::Vec3d(0, 0.f, velocity), particle_generation_rate);
            domain.addFluidSource(src1);
            return;
        case InitialisationType::VortexSmoke:
            velocity = 100;
            concentration = 50;
            return;
        case InitialisationType::CollidingSmoke:
            velocity = 500;
            concentration = 125;
            return;
        case InitialisationType::WaterFall:
            velocity = 0;
            concentration = 70;

            openvdb::BBoxd bbox(i2w_transform->indexToWorld(openvdb::Vec3R(
                                    center_point[0] - 15.5f, center_point[1] - 15.5f, center_point[2] - 15.5f)),
                                i2w_transform->indexToWorld(openvdb::Vec3R(
                                    center_point[0] + 15.5f, center_point[1] + 15.5f, center_point[2] + 15.5f)));

            l = openvdb::tools::createLevelSetBox<openvdb::FloatGrid, openvdb::Vec3d>(bbox, *i2w_transform);
            openvdb::tools::GridSampler<openvdb::FloatGrid::Accessor, openvdb::tools::BoxSampler> sampler(
                l->getAccessor(), l->transform());

            src1 = FluidSource(l->deepCopy(), openvdb::Vec3d(0, 0., -1.f), particle_generation_rate);
            l = openvdb::tools::createLevelSetCube<openvdb::FloatGrid>(0.25f, openvdb::Vec3d(center_point), width / nx);
            l->setTransform(i2w_transform);

            src2 = FluidSource(l->deepCopy(), openvdb::Vec3d(0, 0.f, 0.f), particle_generation_rate);

            domain.addFluidSource(src1);
            // domain.addFluidSource(src2);
            return;
    }
}

int main(int argc, char **argv) {
    std::cout.precision(2);
    // Create a VDB file object.
    openvdb::io::File file("mygrid.vdb");
    // Add the grid pointer to a container.
    openvdb::GridPtrVec grids;

    int frame = 0;
    double voxel_size = width / ny;
    auto i2w_transform = openvdb::math::Transform::createLinearTransform(voxel_size);

    i2w_transform->preTranslate(openvdb::Vec3d(0.5, 0.5, 0.5));

    openvdb::Vec3f center_point(width / 2.f, width / 2.f, width / 2.f);
    auto solidBBox = openvdb::tools::createLevelSetBox<openvdb::FloatGrid, openvdb::Vec3d>(
        openvdb::BBoxd(i2w_transform->indexToWorld(openvdb::Vec3d(1, 1, 1)),
                       i2w_transform->indexToWorld(openvdb::Vec3d(nx - 2, ny - 2, nz - 2))),
        *i2w_transform);
    FluidSimulator sim(solidBBox->evalActiveVoxelBoundingBox());
    FluidDomain fluidDomain(nx, ny, width, width, solidBBox, openvdb::createLevelSet<openvdb::FloatGrid>(voxel_size),
                            i2w_transform, voxel_size);
    fluidDomain.fluidLevelSet().getLevelSet()->setTransform(i2w_transform);
    for (auto iter = fluidDomain.solidLevelSet().getLevelSet()->tree().beginRootChildren(); iter; ++iter) {
        iter->negate();
    }
    fluidDomain.solidLevelSet().getLevelSet()->setTransform(i2w_transform);

    grids.push_back(fluidDomain.solidLevelSet().getLevelSet());
    // for (openvdb::FloatGrid::ValueOnCIter iter = fluidDomain.solidLevelSet()->cbeginValueOn(); iter; ++iter) {
    //     std::cout << "Grid" << iter.getCoord() << " = " << *iter << std::endl;
    // }

    Renderer renderer(nz, ny);
    renderer.initialisation();

    InitialisationType initType = INIT_TYPE;
#if defined EULER_GRID
#elif defined PIC_FLIP
    pic_initialisation(fluidDomain, initType, i2w_transform);
#endif

    const std::filesystem::path vdb_folder{"sequence"};
    if (!std::filesystem::exists(vdb_folder)) std::filesystem::create_directory(vdb_folder);

    sf::RenderWindow window = sf::RenderWindow(sf::VideoMode(WINDOW_SIZE, WINDOW_SIZE), "SFML works!");
    while (window.isOpen() && frame < max_frames_n) {
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
        if (frame % 2000 == 1) {
            fluidDomain.removeFluidSource(0);
            fluidDomain.removeFluidSource(0);
        }

#if defined EULER_GRID
#elif defined PIC_FLIP
        sim.advance_flip_pic(fluidDomain, t_frame, 0.95);
#endif

        grids.push_back(fluidDomain.fluidLevelSet().getLevelSet()->deepCopy());
        std::string filename = "sequence_" + std::to_string(frame) + ".vdb";
        openvdb::io::File((vdb_folder / std::filesystem::path(filename)).string())
            .write({ fluidDomain.fluidLevelSet().getLevelSet()->deepCopy() });

        // grids.push_back(fluidDomain.solidLevelSet().getLevelSet());
        // Sleep(50);
        frame++;
    }
    file.write(grids);
    file.close();
    return 0;
}