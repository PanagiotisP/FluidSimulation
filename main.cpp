// 2D Fluid Simulation for now
// No rectangular bounding box, 4 free surfaces for now
#include "array2.h"
#include "vec.h"
#include "array2_utils.h"
#include "FluidSimulator.h"
#include <SFML/Graphics.hpp>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

#define PI 3.14159265f

// Basic fluid algorithm
//• Start with an initial divergence - free velocity field u0.
//• For time step n = 0, 1, 2, . . .
//• Determine a good time step Dt to go from time t^n to time t^(n + 1).
//• Set u^A = advect(u^n, Dt, u^n).
//• Add u^B = u^A + Dtg.
//• Set u^(n + 1) = project(Dt, u^B).

//flip method
//• Transfer particle values qp to the grid qi, j, k, through equations like(7.1)
//or (7.2), and extrapolate on the grid as necessary.
//• Save the grid values qi, j, k.
//• Compute all other terms on the grid, such as pressure projection, to
//get an updated qnew
//i, j, k.
//• For each particle, interpolate the change qi, j, k = qnew
//i, j, k ? qi, j, k from
//the grid to add it to the particle’s value.
//• Advect the particles in the grid velocity field.

// Time span of each frame
float t_frame = 0.01f;

int max_frames_n = 100;

float density = 1;
// Grid dimensions
int nx = 6;
int ny = nx;
float width = 1;

Array2<Vec2f, Array1<Vec2f>> interpolated_velocity(nx, ny);

void generate_interpolated_velocity(FluidSimulator& sim) {
  for (int j = 0; j < sim.ny; ++j) {
    for (int i = 0; i < sim.nx; ++i) {
      interpolated_velocity(i, j) = sim.evaluate_velocity(Vec2f(i * sim.dx, j * sim.dx));
    }
  }
}

// Axis aligned rectangular Box
//float liquid_bound_x0 = 0, liquid_bound_x1 = width;
//float liquid_bound_y0 = 0, liquid_bound_y1 = width / 2;
//
//float solid_bound_x0 = 0, solid_bound_x1 = width;
//float solid_bound_y0 = 0, solid_bound_y1 = width;


//float distance_from_axis_aligned_rect_box(const float x0, const float x1, const float y0, const float y1, const Vec2f& point) {
//  //inside
//  if (x0 < point[0] && point[0] < x1
//    && y0 < point[1] && point[1] < y1) {
//    return max(x0 - point[0], point[0] - x1, y0 - point[1], point[1] - y1);
//  }
//  //outside
//  else {
//    //closest point to bounding box (p,q)
//    float p, q;
//    if (point[0] < x0) p = x0;
//    else if (x1 < point[0]) p = x1;
//    else p = point[0];
//
//    if (point[1] < y0) q = y0;
//    else if (y1 < point[1]) q = y1;
//    else q = point[1];
//    return sqrt((sqr(point[0] - p)) + (sqr(point[1] - q)));
//  }
//}
//
////Rectangular Box
//float query_liquid_phi(const Vec2f& point) {
//  return distance_from_axis_aligned_rect_box(liquid_bound_x0, liquid_bound_x1, liquid_bound_y0, liquid_bound_y1, point);
//}
//
////Rectangular Box
//float query_boundary_phi(const Vec2f& point) {
//  return -distance_from_axis_aligned_rect_box(solid_bound_x0, solid_bound_x1, solid_bound_y0, solid_bound_y1, point);
//}

//Testing-debugging purpose
void custom_initialisation(FluidSimulator& sim) {
  for (int i = 0; i < sim.nx + 1; ++i) {
    sim.u_grid(i, 0) = 10;
  }
  //sim.u_grid(2, 0) = 10;
  //sim.u_grid(1, 0) = 10;
  for (int j = 0; j < sim.ny + 1; ++j) {
    sim.v_grid(0, j) = 10;
  }
}

int main(int argc, char** argv) {
  int frame = 0;
  float grid_square_width = nx * 10;

  FluidSimulator sim;
  sim.initialize(width, density, nx, ny);
  sim.set_boundary();

  custom_initialisation(sim);

  generate_interpolated_velocity(sim);
  std::cout << "Initial Velocity Field";
  std::cout << interpolated_velocity << std::endl;

  sf::RenderWindow window(sf::VideoMode(800, 800), "SFML works!");
  int grid_width = 500;
  assert(grid_width < 800);


  Array2<sf::RectangleShape> grid(sim.nx, sim.ny);
  Array2<sf::RectangleShape> vectors(sim.nx, sim.ny);

  sf::RectangleShape grid_rectangle(sf::Vector2f(grid_width / sim.nx, grid_width / sim.ny));
  grid_rectangle.setOutlineThickness(-2);
  grid_rectangle.setOutlineColor(sf::Color::Green);

  grid.assign(grid_rectangle);


  sf::RectangleShape vector_arrow(sf::Vector2f(0, 2));
  vector_arrow.setFillColor(sf::Color::Black);

  vectors.assign(vector_arrow);

  sf::Vector2f offset((800 - grid_width) / 2, (800 - grid_width) / 2);

  for (int j = 0; j < sim.ny; ++j)
    for (int i = 0; i < sim.nx; ++i) {
      grid(i, j).setPosition(offset.x + i * grid_width / sim.nx, offset.y + grid_width - (j+1) * grid_width / sim.ny);
      vectors(i, j).setPosition(offset.x + (i + 0.5) * grid_width / sim.nx, offset.y + grid_width - (j + 0.5) * grid_width / sim.ny);
    }

  window.display();
  while (window.isOpen() && (true || frame < 50)) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type == sf::Event::Closed)
        window.close();
    }

    sim.advance(t_frame);
    generate_interpolated_velocity(sim);
    std::cout << "Velocity Field";
    std::cout << interpolated_velocity << std::endl;
    std::cout << "Temperature";
    std::cout<<sim.temperature<<std::endl;

    for (int j = 0; j < sim.ny; ++j) {
      for (int i = 0; i < sim.nx; ++i) {
        Vec2f point(i, j);
        Vec2f velocity_vector = sim.evaluate_velocity(point);
        float angle = atan(fabs(velocity_vector[1]) / fabs(velocity_vector[0])) * 180.f / PI;

        if (velocity_vector[0] >= 0 && velocity_vector[1] >= 0) angle = -angle;
        else if (velocity_vector[0] < 0 && velocity_vector[1] >= 0) angle = -(180 - angle);
        else if (velocity_vector[0] < 0 && velocity_vector[1] < 0) angle = -(180 + angle);
        else if (velocity_vector[0] >= 0 && velocity_vector[1] < 0) angle = angle;

        vectors(i, j).setRotation(angle);
        vectors(i, j).setSize(sf::Vector2f(150 / sim.nx * mag(velocity_vector), 2));
      }
    }

    window.clear(sf::Color::Black);
    for (int i = 0; i < grid.a.size(); ++i) {
      assert(grid.a.size() == vectors.a.size());
      window.draw(grid.a[i]);
      window.draw(vectors.a[i]);
    }

    window.display();
    Sleep(50);
    frame++;
  }

  return 0;
}