#pragma once
#include "FluidSimulator.h"
#include "vec.h"
#include "array2_utils.h"

float liquid_bound_x0, liquid_bound_x1;
float liquid_bound_y0, liquid_bound_y1;

float solid_bound_x0, solid_bound_x1;
float solid_bound_y0, solid_bound_y1;


float distance_from_axis_aligned_rect_box(const float x0, const float x1, const float y0, const float y1, const Vec2f& point) {
  //inside
  if (x0 < point[0] && point[0] < x1
    && y0 < point[1] && point[1] < y1) {
    return max(x0 - point[0], point[0] - x1, y0 - point[1], point[1] - y1);
  }
  //outside
  else {
    //closest point to bounding box (p,q)
    float p, q;
    if (point[0] < x0) p = x0;
    else if (x1 < point[0]) p = x1;
    else p = point[0];

    if (point[1] < y0) q = y0;
    else if (y1 < point[1]) q = y1;
    else q = point[1];
    return sqrt((sqr(point[0] - p)) + (sqr(point[1] - q)));
  }
}

//Rectangular Box
float query_liquid_phi(const Vec2f& point) {
  return distance_from_axis_aligned_rect_box(liquid_bound_x0, liquid_bound_x1, liquid_bound_y0, liquid_bound_y1, point);
}

float query_boundary_phi(const Vec2f& point) {
  return -distance_from_axis_aligned_rect_box(solid_bound_x0, solid_bound_x1, solid_bound_y0, solid_bound_y1, point);
}

FluidSimulator::FluidSimulator() {}

FluidSimulator::~FluidSimulator() {}

void FluidSimulator::initialize(float width, float density, int nx, int ny) {
  this->density = density;
  this->nx = nx;
  this->ny = ny;
  dx = width / nx;

  //TODO Delete these please
  liquid_bound_x0 = 0; liquid_bound_x1 = width - dx / 2;
  liquid_bound_y0 = 0; liquid_bound_y1 = width - dx / 2;

  solid_bound_x0 = 0-dx/2; solid_bound_x1 = width + dx / 2;
  solid_bound_y0 = 0-dx/2; solid_bound_y1 = width + dx / 2;


  u_grid.resize(nx + 1, ny); temp_u_grid.resize(nx + 1, ny);
  v_grid.resize(nx, ny + 1); temp_v_grid.resize(nx, ny + 1);

  temperature.resize(nx, ny); temp_temperature.resize(nx, ny);
  concentration.resize(nx, ny); temp_concentration.resize(nx, ny);

  divirgence.resize(nx, ny);

  solid_phi.resize(nx + 1, ny + 1);
  liquid_phi.resize(nx + 1, ny + 1);

  u_grid.set_zero();
  v_grid.set_zero();
  solid_phi.set_zero();

  temperature.assign(nx, ny, 273 + 27);
}

Vec2f FluidSimulator::evaluate_velocity(const Vec2f& point) {
  float u_value = interpolate_value(point / dx + Vec2f(0.5f, 0), u_grid);
  float v_value = interpolate_value(point / dx + Vec2f(0, 0.5f), v_grid);
  return Vec2f(u_value, v_value);
}

float FluidSimulator::evaluate_scalar(Vec2f& point, Array2f& grid) {
  return interpolate_value(point / dx, grid);
}

void FluidSimulator::set_boundary() {
  for (int j = 0; j < ny + 1; ++j) {
    for (int i = 0; i < nx + 1; ++i) {
      solid_phi(i, j) = query_boundary_phi(Vec2f((i+0.5) * dx, (j+0.5) * dx));
      liquid_phi(i, j) = query_liquid_phi(Vec2f((i+0.5) * dx, (j+0.5) * dx));
    }
  }
  //std::cout << solid_phi << std::endl;
  Array2f interpolated_liquid_phi(nx, ny);
  Vec2f offset(0.5, 0.5);
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      interpolated_liquid_phi(i,j) = interpolate_value(offset+Vec2f(i, j), liquid_phi);
    }
  }
  std::cout << interpolated_liquid_phi << std::endl;
}

void FluidSimulator::advance(float t_frame) {
  float t = 0;
  while (t < t_frame) {
    float timestep = compute_cfl();
    //TODO DELETE ME TEST VALUE
    timestep = 0.01f;
    if (timestep + t >= t_frame) {
      timestep = t_frame - t;
    }
    advect(timestep);
    //apply forces (e.g. gravity)
    add_forces(timestep);

    Array2f divirgence(nx, ny);
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {

      }
    }


    //update temperature and concentration
    float rate_t = 10;
    float rate_s = 1;
    float target_temperature = 273 + 37;
    //temperature(nx / 2, ny / 2) += (1 - exp(-rate_t * timestep)) * (target_temperature - temperature(nx / 2, ny / 2));
    temperature(0, 0) += (1 - exp(-rate_t * timestep)) * (target_temperature - temperature(0, 0));
    temperature(0, 1) += (1 - exp(-rate_t * timestep)) * (target_temperature - temperature(0, 1));
    concentration(nx / 2, ny / 2) += clamp(rate_s * timestep, 0.f, 1.f);

    project(timestep);
    t += timestep;
  }
}

Vec2f FluidSimulator::trace_rk2(const Vec2f& point, float dt) {
  Vec2f mid_pos = point - (float)0.5 * dt * evaluate_velocity(point);
  return point - dt * evaluate_velocity(mid_pos);
}

float FluidSimulator::compute_cfl() {
  float max_vel = 0;

  for (int i = 0; i < u_grid.a.size(); ++i) {
    max_vel = max(u_grid.a[i], fabs(max_vel));
  }

  for (int i = 0; i < u_grid.a.size(); ++i) {
    max_vel = max(u_grid.a[i], fabs(max_vel));
  }

  for (int i = 0; i < u_grid.a.size(); ++i) {
    max_vel = max(u_grid.a[i], fabs(max_vel));
  }
  max_vel += sqrt(5 * dx * grav);
  return dx / max_vel;
}

void FluidSimulator::advect(float dt) {
  temp_u_grid.set_zero();
  temp_v_grid.set_zero();
  temp_temperature.set_zero();
  temp_concentration.set_zero();

  //hypothetical particle for u-component
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx + 1; ++i) {
      Vec2f pos((i - 0.5) * dx, j * dx);
      pos = trace_rk2(pos, dt);
      temp_u_grid(i, j) = evaluate_velocity(pos)[0];
    }
  }

  //hypothetical particle for v-component
  for (int j = 0; j < ny + 1; ++j) {
    for (int i = 0; i < nx; ++i) {
      Vec2f pos(i * dx, (j - 0.5) * dx);
      pos = trace_rk2(pos, dt);
      temp_v_grid(i, j) = evaluate_velocity(pos)[1];
    }
  }

  // Advect other quantities like Temperature and Concentration
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      Vec2f pos(i * dx, j * dx);
      Vec2f starting_pos = trace_rk2(pos, dt);
      temp_temperature(i, j) = evaluate_scalar(starting_pos, temperature);
      temp_concentration(i, j) = evaluate_scalar(starting_pos, concentration);
    }
  }

  temperature = temp_temperature;
  concentration = temp_concentration;
  u_grid = temp_u_grid;
  v_grid = temp_v_grid;

  calculate_divirgence();
  std::cout << "Divirgence";
  std::cout << divirgence << std::endl;

}

void FluidSimulator::add_forces(float dt) {
  //gravity
  for (int j = 0; j < ny + 1; ++j)
    for (int i = 0; i < nx; ++i) {
      // Temporarily moving the grid perpendicular to the vertical axis
      v_grid(i, j) -= 0;//grav*dt;
    }
}

void FluidSimulator::project(float dt) {
  const int system_size = nx * ny;

  if (rhs.size() != system_size) {
    rhs.resize(system_size);
    pressure_grid.resize(system_size);
    alpha_matrix.resize(system_size);
  }

  alpha_matrix.zero();
  rhs.assign(rhs.size(), 0);
  pressure_grid.assign(pressure_grid.size(), 0);


  float scale = dt / (density * sqr(dx));
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int idx = i + j * nx;

      rhs[idx] = 0;
      //pressure_grid[idx] = 0;
      //solid velocity zero for now
      if (query_liquid_phi(Vec2f(i * dx, j * dx)) <= 0) {
        //left neighbour
        if (query_liquid_phi(Vec2f((i - 1) * dx, j * dx)) <= 0) {
          alpha_matrix.add_to_element(idx, idx, scale);
          alpha_matrix.add_to_element(idx, idx - 1, -scale);
        }
        else
          alpha_matrix.add_to_element(idx, idx, scale);

        rhs[idx] += u_grid(i, j) / dx;
        if (query_boundary_phi(Vec2f(i - 1, j)) <= 0) {
          rhs[idx] -= u_grid(i, j) / dx;
        }

        //right neighbour
        if (query_liquid_phi(Vec2f((i + 1) * dx, j * dx)) <= 0) {
          alpha_matrix.add_to_element(idx, idx, scale);
          alpha_matrix.add_to_element(idx, idx + 1, -scale);
        }
        else
          alpha_matrix.add_to_element(idx, idx, scale);
        rhs[idx] -= u_grid(i + 1, j) / dx;
        if (query_boundary_phi(Vec2f(i + 1, j)) <= 0) {
          rhs[idx] += u_grid(i + 1, j) / dx;
        }

        //bottom neighbour
        if (query_liquid_phi(Vec2f(i * dx, (j - 1) * dx)) <= 0) {
          alpha_matrix.add_to_element(idx, idx, scale);
          alpha_matrix.add_to_element(idx, idx - nx, -scale);
        }
        else
          alpha_matrix.add_to_element(idx, idx, scale);
        rhs[idx] += v_grid(i, j) / dx;
        if (query_boundary_phi(Vec2f(i, j - 1)) <= 0) {
          rhs[idx] -= v_grid(i, j) / dx;
        }
        //top neighbour
        if (query_liquid_phi(Vec2f(i * dx, (j + 1) * dx)) <= 0) {
          alpha_matrix.add_to_element(idx, idx, scale);
          alpha_matrix.add_to_element(idx, idx + nx, -scale);
        }
        else
          alpha_matrix.add_to_element(idx, idx, scale);
        rhs[idx] -= v_grid(i, j + 1) / dx;
        if (query_boundary_phi(Vec2f(i, j + 1)) <= 0) {
          rhs[idx] += v_grid(i, j + 1) / dx;
        }
      }
    }
  }
  //Solve the system using Robert Bridson's incomplete Cholesky PCG solver

  double tolerance;
  int iterations;
  solver.set_solver_parameters(1e-18, 1000);
  bool success = solver.solve(alpha_matrix, rhs, pressure_grid, tolerance, iterations);
  printf("Solver took %d iterations and had residual %e\n", iterations, tolerance);
  if (!success) {
    printf("WARNING: Pressure solve failed!************************************************\n");
  }

  for (int j = 0; j < ny; ++j) {
    for (int i = 1; i < nx; ++i) {
      int idx = i + nx * j;
      u_grid(i, j) -= dt * (float)(pressure_grid[idx] - pressure_grid[idx-1]) / dx;
    }
  }

  for (int j = 1; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      int idx = i + nx * j;
      v_grid(i, j) -= dt * (float)(pressure_grid[idx] - pressure_grid[idx-nx]) / dx;
    }
  }
}

void FluidSimulator::calculate_divirgence() {
  float scale = 1 / dx;
  for (int j = 0; j < ny; ++j) {
    for (int i = 0; i < nx; ++i) {
      divirgence(i, j) = scale * (u_grid(i + 1, j) - u_grid(i, j)
        + v_grid(i, j + 1) - v_grid(i, j));
    }
  }
}
