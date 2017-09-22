#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <utility>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//

// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// Reference Stuff
const double ref_cte = 0;
const double ref_epsi = 0;
const double ref_v = 40;

// Costs
const double cte_cost = 20;
const double epsi_cost = 15;
const double velocity_cost = 1;
const double delta_cost = 5;
const double a_cost = 1;
const double delta_d_cost = 800;
const double delta_a_cost = 5;

// Indices
const size_t x_start = 0;
const size_t y_start = x_start + N;
const size_t psi_start = y_start + N;
const size_t v_start = psi_start + N;
const size_t cte_start = v_start + N;
const size_t epsi_start = cte_start + N;
const size_t delta_start = epsi_start + N;
const size_t a_start = delta_start + N - 1;


class FG_eval {
public:
    // Fitted polynomial coefficients
    Eigen::VectorXd coeffs;

    explicit FG_eval(Eigen::VectorXd coefficients) :
        coeffs(std::move(coefficients)) {}

    using ADvector = CPPAD_TESTVECTOR(AD<double>);

    void operator()(ADvector &fg, const ADvector &vars) {
      // TODO: implement MPC
      // `fg` a vector of the cost constraints, `vars` is a vector of variable
      // values (state & actuators)
      // NOTE: You'll probably go back and forth between this function and
      // the Solver function below.

      // The cost is stored is the first element of `fg`.
      // Any additions to the cost should be added to `fg[0]`.
      fg[0] = 0;

      // The part of the cost based on the reference state.
      for (size_t t = 0; t < N; t++) {
        fg[0] += cte_cost * CppAD::pow(vars[cte_start + t] - ref_cte, 2);
        fg[0] += epsi_cost * CppAD::pow(vars[epsi_start + t] - ref_epsi, 2);
        fg[0] += velocity_cost * CppAD::pow(vars[v_start + t] - ref_v, 2);
      }

      // Minimize the use of actuators.
      for (size_t t = 0; t < N - 1; t++) {
        fg[0] += delta_cost * CppAD::pow(vars[delta_start + t], 2);
        fg[0] += a_cost * CppAD::pow(vars[a_start + t], 2);
      }

      // Minimize the value gap between sequential actuations.
      for (size_t t = 0; t < N - 2; t++) {
        fg[0] += delta_d_cost *
                 CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t],
                            2);
        fg[0] += delta_a_cost *
                 CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
      }

      // Initial constraints
      //
      // We add 1 to each of the starting indices due to cost being located at
      // index 0 of `fg`.
      // This bumps up the position of all the other values.

      fg[1 + x_start] = vars[x_start];
      fg[1 + y_start] = vars[y_start];
      fg[1 + psi_start] = vars[psi_start];
      fg[1 + v_start] = vars[v_start];
      fg[1 + cte_start] = vars[cte_start];
      fg[1 + epsi_start] = vars[epsi_start];

      // The rest of the constraints
      for (size_t t1 = 1; t1 < N; t1++) {
        const size_t t0 = t1 - 1;

        // The state at time t.
        const AD<double> x0 = vars[x_start + t0];
        const AD<double> y0 = vars[y_start + t0];
        const AD<double> psi0 = vars[psi_start + t0];
        const AD<double> v0 = vars[v_start + t0];
        const AD<double> cte0 = vars[cte_start + t0];
        const AD<double> epsi0 = vars[epsi_start + t0];

        // The state at time t+1 .
        const AD<double> x1 = vars[x_start + t1];
        const AD<double> y1 = vars[y_start + t1];
        const AD<double> psi1 = vars[psi_start + t1];
        const AD<double> v1 = vars[v_start + t1];
        const AD<double> cte1 = vars[cte_start + t1];
        const AD<double> epsi1 = vars[epsi_start + t1];

        // Only consider the actuation at time t.
        const AD<double> delta0 = vars[delta_start + t0];
        const AD<double> a0 = vars[a_start + t0];
        const AD<double> x0_2 = x0 * x0;
        const AD<double> x0_3 = x0 * x0 * x0;
        const AD<double> f0 =
            coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0_2 + coeffs[3] * x0_3;
        const AD<double> psides0 =
            CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * x0_2);
        fg[1 + x_start + t1] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
        fg[1 + y_start + t1] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
        fg[1 + psi_start + t1] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
        fg[1 + v_start + t1] = v1 - (v0 + a0 * dt);
        fg[1 + cte_start + t1] =
            cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
        fg[1 + epsi_start + t1] =
            epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
      }

    }
};

//
// MPC class definition implementation.
//
MPC::MPC() = default;

MPC::~MPC() = default;

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  const size_t n_state_vars = 6;
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];
  double steering_angle = state[6];
  double throttle = state[7];
  const size_t n_actuator_vars = 2;
  size_t n_vars = n_state_vars * N + n_actuator_vars * (N - 1);

  // TODO: Set the number of constraints
  size_t n_constraints = N * 6;

  // Initial value of the independent variables.
  // Should be 0 except for the initial values.
  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0.0;
  }
  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (size_t i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }
  int delay_indices = 1;
  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (size_t i = delta_start; i < delta_start + delay_indices; i++) {
    // Initialize with current value.
    // Because we get stuff with 100ms delay.
    vars_lowerbound[i] = steering_angle;
    vars_upperbound[i] = steering_angle;
  }

  for (size_t i = delta_start + delay_indices; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332 * Lf;
    vars_upperbound[i] = 0.436332 * Lf;
  }

  // Acceleration/decceleration upper and lower limits.
  for (size_t i = a_start; i < a_start + delay_indices; i++) {
    // Initialize with current value.
    // Because we get stuff with 100ms delay.
    vars_lowerbound[i] = throttle;
    vars_upperbound[i] = throttle;
  }
  for (size_t i = a_start + delay_indices; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for constraints
  // All of these should be 0 except the initial
  // state indices.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (size_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(std::move(coeffs));

  // options for IPOPT solver
  std::string options;
  options += "Integer print_level  0\n";
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  options += "Numeric max_cpu_time          0.5\n";

  // Place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  bool ok = solution.status == CppAD::ipopt::solve_result<Dvector>::success;
  if (!ok) {
    throw "solution not found";
  }

  // Cost
  std::vector<double> result;
  result.push_back(solution.x[delta_start + delay_indices]);
  result.push_back(solution.x[a_start + delay_indices]);
  for (size_t i = 1; i < N; ++i) {
    result.push_back(solution.x[x_start + i]);
    result.push_back(solution.x[y_start + i]);
  }
  return result;
}
