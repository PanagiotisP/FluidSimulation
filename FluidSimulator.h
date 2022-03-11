#pragma once
#include "array2.h"
#include "pcgsolver/pcg_solver.h"
#include "pcgsolver/sparse_matrix.h"
#include "vec.h"
class FluidSimulator
{
public:
  void initialize(float width, float density, int nx, int ny);
  //void set_boundary(float (*phi)(const Vec3f&));
  //void set_liquid(float (*phi)(const Vec3f&));
  //void add_particle(const Vec3f& pos);
  
  //Get velocity at point
  //Coordinates given on original grid
  Vec2f evaluate_velocity(const Vec2f& point);
  
  //Get scalar quantity at point
  float evaluate_scalar(Vec2f& point, Array2f& grid);

  void set_boundary();

  void advance(float t_frame);

  //Grid dimensions
  int nx, ny;

  float dx;

  const float grav = 9.81f;
  float density;

  //Boundaries
  Array2f solid_phi;
  Array2f liquid_phi;

  //Velocity vector fields
  Array2f u_grid, temp_u_grid;
  Array2f v_grid, temp_v_grid;

  //Scalar fields for other quantities
  Array2f temperature, temp_temperature;
  Array2f concentration, temp_concentration;

  Array2f divirgence;

  //Pressure solver
  std::vector<double> rhs;
  std::vector<double> pressure_grid;
  SparseMatrixd alpha_matrix;
  PCGSolver<double> solver;

  FluidSimulator();
  ~FluidSimulator();

private:
  Vec2f trace_rk2(const Vec2f& point, float dt);

  float compute_cfl();

  //void advect_particles(float dt);
  //Semi-Lagrangian advection
  void advect(float dt);

  //Add external forces like gravity
  void add_forces(float dt);

  //Apply incrompressibility constraint
  void project(float dt);
  //void constrain_velocity();

  //helpers for pressure projection
  //void compute_weights();
  //void solve_pressure(float dt);
  //void compute_phi

  void calculate_divirgence();
};
