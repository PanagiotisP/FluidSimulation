// 2D Fluid Simulation for now
// No rectangular bounding box, 4 free surfaces for now

#include "array2.h"
#include "vec.h"
#include "array2_utils.h"
#include "pcgsolver/sparse_matrix.h"
#include "pcgsolver/pcg_solver.h"
#include <SFML/Graphics.hpp>
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

#define PI 3.14159265

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
bool end_of_frame = false;


float grav = 9.81f;
float density = 1;
// Grid dimensions
int nx = 2;
int ny = 2;
float width = nx;
float dx = width/float(nx);

Array2f u_grid(nx+1, ny), temp_u_grid(nx+1, ny);
Array2f v_grid(nx, ny+1), temp_v_grid(nx, ny+1);
Array2<Vec2f, Array1<Vec2f>> interpolated_velocity(nx, ny);
Array2f temperature(nx, ny), temp_temperature(nx, ny);
Array2f concentration(nx, ny), temp_concentration(nx, ny);

std::vector<double> rhs;
std::vector<double> pressure_grid;
SparseMatrixd alpha_matrix;
PCGSolver<double> solver;


Vec2f evaluate_velocity(Vec2f point);

//Rectangular Box
float query_liquid_phi() {
  
}
//
//float query_boundary_phi() {
//  ;
//}

Vec2f runge_kutta(Vec2f starting_pos, float timestep) {
  Vec2f mid_pos = starting_pos-(float)0.5*timestep*evaluate_velocity(starting_pos);
  return starting_pos-timestep*evaluate_velocity(mid_pos);
}

Vec2f evaluate_velocity(Vec2f point) {
  float u_value = interpolate_value(point/dx+Vec2f(0.5f, 0), u_grid);
  float v_value = interpolate_value(point/dx+Vec2f(0, 0.5f), v_grid);
  return Vec2f(u_value, v_value);
}

float evaluate_scalar(Vec2f point, Array2f& grid) {
  return interpolate_value(point/dx, grid);
}

float compute_cfl() {
  float max_vel = 0;

  for (int i = 0; i<u_grid.a.size(); ++i) {
    max_vel = max(u_grid.a[i], fabs(max_vel));
  }

  for (int i = 0; i<u_grid.a.size(); ++i) {
    max_vel = max(u_grid.a[i], fabs(max_vel));
  }

  for (int i = 0; i<u_grid.a.size(); ++i) {
    max_vel = max(u_grid.a[i], fabs(max_vel));
  }
  max_vel += sqrt(5*dx*grav);
  return dx/max_vel;
}

//Semi-Lagrangian advection
void advect(float timestep) {
  temp_u_grid.set_zero();
  temp_v_grid.set_zero();
  temp_temperature.set_zero();
  temp_concentration.set_zero();

  //hypothetical particle for u-component
  for (int i = 0; i<nx+1; ++i) {
    for (int j = 0; j<ny; ++j) {
      Vec2f pos(i*dx, j*dx);
      pos = runge_kutta(pos, timestep);
      temp_u_grid(i, j) = evaluate_velocity(pos)[0];
    }
  }

  //hypothetical particle for v-component
  for (int i = 0; i<nx; ++i) {
    for (int j = 0; j<ny+1; ++j) {
      Vec2f pos(i*dx, j*dx);
      pos = runge_kutta(pos, timestep);
      temp_v_grid(i, j) = evaluate_velocity(pos)[1];
    }
  }

  // Advect other quantities like Temperature and Concentration
  for (int i = 0; i<nx; ++i) {
    for (int j = 0; j<ny; ++j) {
      Vec2f pos(i*dx, j*dx);
      Vec2f starting_pos = runge_kutta(pos, timestep);
      temp_temperature(i, j) = evaluate_scalar(starting_pos, temperature);
      temp_concentration(i, j) = evaluate_scalar(starting_pos, concentration);
    }
  }

  temperature = temp_temperature;
  concentration = temp_concentration;
  u_grid = temp_u_grid;
  v_grid = temp_v_grid;
}

void add_force(float timestep) {
  //gravity
  for (int j = 0; j<ny+1; ++j)
    for (int i = 0; i<nx; ++i) {
      // Temporarily moving the grid perpendicular to the vertical axis
      v_grid(i, j) -= 0;//grav*timestep;
    }
}

void project(float timestep) {
  int system_size = nx*ny;

  if (rhs.size()!=system_size) {
    rhs.resize(system_size);
    pressure_grid.resize(system_size);
    alpha_matrix.resize(system_size);
  }

  alpha_matrix.zero();
  rhs.assign(rhs.size(), 0);
  pressure_grid.assign(pressure_grid.size(), 0);


  float scale = timestep/(density*sqr(dx));
  for (int j = 1; j<ny-1; ++j) {
    for (int i = 1; i<nx-1; ++i) {
      int idx = i+j*nx;

      rhs[idx] = 0;
      pressure_grid[idx] = 0;
      //if cell i,j fluid
      {
        //left neighbour
        //if cell fluid
        alpha_matrix.add_to_element(idx, idx, scale);
        alpha_matrix.add_to_element(idx, idx-1, -scale);
        //else if empty
        //alpha_matrix.add_to_element(idx, idx, scale);
        rhs[idx] += u_grid(i, j)/dx;


        //right neighbour
        //if cell fluid
        alpha_matrix.add_to_element(idx, idx, scale);
        alpha_matrix.add_to_element(idx, idx+1, -scale);
        //else if empty
        //alpha_matrix.add_to_element(idx, idx, scale);
        rhs[idx] -= u_grid(i+1, j)/dx;

        //bottom neighbour
        //if cell fluid
        alpha_matrix.add_to_element(idx, idx, scale);
        alpha_matrix.add_to_element(idx, idx-nx, -scale);
        //else if empty
        //alpha_matrix.add_to_element(idx, idx, scale);
        rhs[idx] += u_grid(i, j)/dx;

        //top neighbour
        //if cell fluid
        alpha_matrix.add_to_element(idx, idx, scale);
        alpha_matrix.add_to_element(idx, idx+nx, -scale);
        //else if empty
        //alpha_matrix.add_to_element(idx, idx, scale);
        rhs[idx] -= u_grid(i, j+1)/dx;
      }
    }
  }
  //Solve the system using Robert Bridson's incomplete Cholesky PCG solver

  double tolerance;
  int iterations;
  solver.set_solver_parameters(1e-18, 1000);
  bool success = solver.solve(alpha_matrix, rhs, pressure_grid, tolerance, iterations);
  //printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
  if (!success) {
    printf("WARNING: Pressure solve failed!************************************************\n");
  }

  for (int j = 0; j<ny; ++j) {
    for (int i = 1; i<nx-1; ++i) {
      int idx = i+nx*j;
      u_grid(i, j) -= timestep*(float)(pressure_grid[idx]-pressure_grid[idx-1])/dx;
    }
  }

  for (int j = 1; j<ny-1; ++j) {
    for (int i = 0; i<nx; ++i) {
      int idx = i+nx*j;
      v_grid(i, j) -= timestep*(float)(pressure_grid[idx]-pressure_grid[idx-nx])/dx;
    }
  }
}

void advance(float t_frame) {
  float t = 0;
  while (t<t_frame) {
    float timestep = compute_cfl();
    //TODO DELETE ME TEST VALUE
    timestep = 0.01f;
    if (timestep+t>=t_frame) {
      timestep = t_frame-t;
    }
    advect(timestep);
    //apply forces (gravity)
    add_force(timestep);

    //update temperature and concentration
    float rate_t = 5;
    float rate_s = 1;
    float target_temperature = 273+30;
    temperature(nx/2, ny/2) += (1-exp(-rate_t*timestep))*(target_temperature-temperature(nx/2, ny/2));
    concentration(nx/2, ny/2) += clamp(rate_s*timestep, 0.f, 1.f);

    project(timestep);
    
    for (int j = 0; j<ny; ++j) {
      for (int i = 0; i<nx; ++i) {
        interpolated_velocity(i, j) = evaluate_velocity(Vec2f(i*dx, j*dx));
      }
    }

    t += timestep;
  }
}

void initialize() {
  u_grid.set_zero();
  v_grid.set_zero();

  for (int i = 0; i<interpolated_velocity.a.size(); ++i) {
    interpolated_velocity.a[i] = Vec2f(0, 0);
  }

  for (int j = 0; j<ny; ++j) {
    u_grid(0, j) = 1;
    u_grid(1, j) = 0;
    u_grid(2, j) = 0;
  }
  for(int i =0 ; i<nx;++i){
    v_grid(i, 0) = 0;
    v_grid(i, 1) = 1;
    v_grid(i, 2) = 0;
  }
  //for (int i = 0; i<ny+1; ++i) {
  //  v_grid(0, i) = 10;
  //}

  temperature.assign(nx, ny, 273+15);
}

//void simulate(Array2f u_0, Array2f v_0, Array2f w_0) {
//  int frame_n = 0;
//  while (frame_n< max_frames_n) {
//    select_dt();
//  }
//}

int main(int argc, char** argv) {
  int frame=0;
  float grid_square_width = nx*10;

  initialize();


  for (int j = 0; j<ny; ++j) {
    for (int i = 0; i<nx; ++i) {
      interpolated_velocity(i, j) = evaluate_velocity(Vec2f(i*dx, j*dx));
    }
  }
  std::cout<<"Initial Velocity Field";
  std::cout<<interpolated_velocity<<std::endl;

  sf::RenderWindow window(sf::VideoMode(800, 800), "SFML works!");
  int grid_width = 500;
  assert(grid_width<800);


  Array2<sf::RectangleShape> grid(nx, ny);
  Array2<sf::RectangleShape> vectors(nx, ny);

  sf::RectangleShape grid_rectangle(sf::Vector2f(grid_width/nx, grid_width/ny));
  grid_rectangle.setOutlineThickness(-2);
  grid_rectangle.setOutlineColor(sf::Color::Green);

  grid.assign(grid_rectangle);


  sf::RectangleShape vector_arrow(sf::Vector2f(0, 2));
  vector_arrow.setFillColor(sf::Color::Black);

  vectors.assign(vector_arrow);

  sf::Vector2f offset((800-grid_width)/2, (800-grid_width)/2);

  for (int j = 0; j<ny; ++j)
    for (int i = 0; i<nx; ++i) {
      grid(i, j).setPosition(offset.x+i*grid_width/nx, offset.y+j*grid_width/ny);
      vectors(i, j).setPosition(offset.x+(i+0.5)*grid_width/nx, offset.y+(j+0.5)*grid_width/ny);
    }

  window.display();
  while (window.isOpen()&&(true || frame<50)) {
    sf::Event event;
    while (window.pollEvent(event)) {
      if (event.type==sf::Event::Closed)
        window.close();
    }

    advance(t_frame);
    std::cout<<"Velocity Field";
    std::cout<<interpolated_velocity<<std::endl;
    //std::cout<<temperature<<std::endl;

    for (int j = 0; j<ny; ++j) {
      for (int i = 0; i<nx; ++i) {
        Vec2f velocity_vector = interpolated_velocity(i, j);
        float angle = atan(fabs(velocity_vector[1])/fabs(velocity_vector[0]))*180/PI;
        
        if (velocity_vector[0]>=0&&velocity_vector[1]>=0) angle = -angle;
        else if (velocity_vector[0]<0&&velocity_vector[1]>=0) angle = -(180-angle);
        else if (velocity_vector[0]<0&&velocity_vector[1]<0) angle = -(180+angle);
        else if (velocity_vector[0]>=0&&velocity_vector[1]<0) angle = angle;

        vectors(i, j).setRotation(angle);
        vectors(i, j).setSize(sf::Vector2f(150/nx*mag(velocity_vector), 2));
      }
    }

    window.clear(sf::Color::Black);
    for (int i = 0; i<grid.a.size(); ++i) {
      assert(grid.a.size()==vectors.a.size());
      window.draw(grid.a[i]);
      window.draw(vectors.a[i]);
    }

    window.display();
    //Sleep(1000);
    frame++;
  }

  return 0;
}