#include "FluidDomain.h"
#include "FluidSimulator.h"
#include "FluidSource.h"
#include "Renderer.h"
#include "render_utils.h"
#include "vec.h"

#include <SFML/Graphics.hpp>
#include <array>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <openvdb/tools/FastSweeping.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/LevelSetAdvect.h>
#include <openvdb/tools/LevelSetPlatonic.h>
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/LevelSetUtil.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/ParticlesToLevelSet.h>

#ifdef _WIN32
    #include <Windows.h>
#else
    #include <unistd.h>
#endif

// #define EULER_GRID
#ifndef EULER_GRID
    #define PIC_FLIP
#endif

#define INIT_TYPE InitialisationType::Test

// Time span of each frame
float t_frame = 0.01f;

// Number of frames of the simulation
int max_frames_n = 500;

// Grid dimensions
int nx = 64;
int ny = 64;
int nz = 64;
float width = 1;

enum class InitialisationType { CenterSmoke, VortexSmoke, CollidingSmoke, WaterFall, Test };

// This functions initialises the domain with fluids sources and solid objects
void domain_initialisation(FluidDomain &domain, InitialisationType type, openvdb::math::Transform::Ptr i2w_transform,
                           float voxelSize) {
    MacGrid &grid = domain.grid();
    openvdb::Vec3f center_point(width / 2, width / 2, width / 2);

    // Solid object
    // auto solidBBox = openvdb::tools::createLevelSetCube<openvdb::FloatGrid>(width, center_point, voxel_size, 5);
    auto solidObject1(std::make_shared<SolidObject>(0.1f, voxelSize, 250., i2w_transform, center_point,
                                                    openvdb::Vec3d(0, 0, 0), openvdb::Vec3d(0)));
    auto solidObject2(std::make_shared<SolidObject>(0.1f, voxelSize, 250., i2w_transform,
                                                    center_point + openvdb::Vec3d(0, 0 * width / 4, 0.1),
                                                    openvdb::Vec3d(0, 0, 0), openvdb::Vec3d(0)));
    domain.addSolidObj(solidObject1);
    // domain.addSolidObj(solidObject2);

    // Fluid Sources
    // particles per grid cell per second
    int particle_generation_rate = 625;

    std::shared_ptr<FluidSource> src1;
    std::shared_ptr<FluidSource> src2;
    std::shared_ptr<FluidSource> src3;
    openvdb::FloatGrid::Ptr l;

    switch (type) {
        case InitialisationType::CenterSmoke: break;
        case InitialisationType::VortexSmoke: break;
        case InitialisationType::CollidingSmoke: break;
        case InitialisationType::WaterFall: break;
        case InitialisationType::Test:
            l = openvdb::tools::createLevelSetCube<openvdb::FloatGrid>(
                0.5f, center_point - openvdb::Vec3d(0, width / 4, 0), voxelSize);
            l->setTransform(i2w_transform);
            src1 = std::make_shared<FluidSource>(l->deepCopy(), openvdb::Vec3d(0, 0., 0.f), particle_generation_rate,
                                                 domain.solidLevelSet(), t_frame);

            l = openvdb::tools::createLevelSetCube<openvdb::FloatGrid>(
                0.5f, center_point - openvdb::Vec3d(0, width / 4, -3 * width / 4), voxelSize);
            l->setTransform(i2w_transform);
            src2 = std::make_shared<FluidSource>(l->deepCopy(), openvdb::Vec3d(0, 0.f, 0.f), particle_generation_rate,
                                                 domain.solidLevelSet(), t_frame);

            domain.addFluidSource(src1);
            // domain.addFluidSource(src2);
            return;
    }
}

int main(int argc, char **argv) {
    std::cout.precision(4);
    openvdb::initialize();

    int frame = 1;
    double voxel_size = width / max(nx, ny, nz);
    auto i2w_transform = openvdb::math::Transform::createLinearTransform(voxel_size);
    i2w_transform->preTranslate(openvdb::Vec3d(0.5, 0.5, 0.5));
    openvdb::Vec3f center_point(width / 2.f);

    // Create a VDB file object and grid pointers for visualisation
    openvdb::io::File file1("fluid.vdb");
    openvdb::io::File file2("solid.vdb");
    openvdb::io::File file3("combined.vdb");

    openvdb::GridPtrVec fluid_level_set;
    openvdb::GridPtrVec solid_level_set;
    openvdb::GridPtrVec combined_level_set;

    // Solid boundary
    // auto solidBBox = openvdb::tools::createLevelSetCube<openvdb::FloatGrid>(width, center_point, voxel_size, 5);
    auto solidBBox = openvdb::tools::createLevelSetBox<openvdb::FloatGrid, openvdb::Vec3d>(
        openvdb::BBoxd(i2w_transform->indexToWorld(openvdb::Vec3d(nx / 4, 1, nz / 4)),
                       i2w_transform->indexToWorld(openvdb::Vec3d(3 * nx / 4, ny - 2, 3 * nz / 4))),
        *i2w_transform, 5);
    solidBBox->setTransform(i2w_transform);

    FluidSimulator sim(openvdb::CoordBBox(openvdb::Coord(0, 0, 0), openvdb::Coord(nx, ny, nz)));

    FluidDomain fluidDomain(solidBBox, i2w_transform, voxel_size);

    // Negate solidBBox's level set to turn it inside-out
    fluidDomain.solidLevelSet().invert();

    fluid_level_set.push_back(fluidDomain.solidLevelSet().getLevelSet());

    Renderer renderer(nz, ny);
    renderer.initialisation();

    InitialisationType initType = INIT_TYPE;
#if defined EULER_GRID
#elif defined PIC_FLIP
    domain_initialisation(fluidDomain, initType, i2w_transform, voxel_size);
#endif

    const std::filesystem::path vdb_folder { "sequence" };
    if (!std::filesystem::exists(vdb_folder)) std::filesystem::create_directory(vdb_folder);

    // Renderer window
    auto window = sf::RenderWindow(sf::VideoMode(WINDOW_SIZE, WINDOW_SIZE), "SFML works!");
    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = t_start;

    while (window.isOpen() && frame <= max_frames_n) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) window.close();
        }

        window.clear(sf::Color::Black);
        renderer.draw(fluidDomain, window);
        window.display();

        auto t_frame_start = std::chrono::high_resolution_clock::now();
#if defined EULER_GRID
#elif defined PIC_FLIP
        sim.advance_flip_pic(fluidDomain, t_frame, 0.95);
#endif
        t_end = std::chrono::high_resolution_clock::now();
        printf("\nFrame %d. Last frame time %.3f. Average time per frame %.3f\n", frame,
               std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_frame_start).count() / 1000.0,
               std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() / 1000.0 / frame);

        // Store current frame's level sets to vector
        fluid_level_set.push_back(fluidDomain.fluidLevelSet().getLevelSet()->deepCopy());
        if (fluidDomain.solidObjs().size() > 0) {
            solid_level_set.push_back(fluidDomain.getSolidObj(0)->levelSet().getLevelSet()->deepCopy());
            combined_level_set.push_back(
                openvdb::tools::csgUnionCopy(*fluidDomain.fluidLevelSet().getLevelSet()->deepCopy(),
                                             *fluidDomain.getSolidObj(0)->levelSet().getLevelSet()->deepCopy()));
        }

        // std::string filename = "sequence_" + std::to_string(frame) + ".vdb";
        // openvdb::io::File((vdb_folder / std::filesystem::path(filename)).string())
        //     .write({ fluidDomain.fluidLevelSet().getLevelSet()->deepCopy() });

        // fluid_level_set.push_back(fluidDomain.solidLevelSet().getLevelSet());
        // Sleep(50);
        frame++;
    }

    // Write and close files
    file1.write(fluid_level_set);
    file2.write(solid_level_set);
    file3.write(combined_level_set);

    file1.close();
    file2.close();
    file3.close();
    return 0;
}