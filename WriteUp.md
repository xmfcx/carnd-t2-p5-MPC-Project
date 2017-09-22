# Implementation

> Student describes their model in detail. This includes the state, actuators and update equations.

In this model, we are trying to minimize the cost function that is determined by 6 state variables:

* x : x position
* y : y position
* psi : orientation
* v : velocity
* cte : cross track error
* epsi : orientation error

Also there are 2 actuator variables:

* a : acceleration
* delta : steering angle

Here delta is limited between -25 and +25 degrees.

Also the cost of cte and epsi is larger than the cost of velocity. 
Because we value staying on track higher than speeding.
Also in the same way, higher delta_cost compared to acceleration cost.

And the model I used is bicycle model from the Putting it all together lesson: https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/f1820894-8322-4bb3-81aa-b26b3c6dcbaf/lessons/338b458f-7ebf-449c-9ad1-611eb933b076/concepts/d26b8460-653f-4479-bc24-68bb62c146ba

And update functions are modified from the solution of the Mind the line lesson: https://classroom.udacity.com/nanodegrees/nd013/parts/40f38239-66b6-46ec-ae68-03afd8a601c8/modules/f1820894-8322-4bb3-81aa-b26b3c6dcbaf/lessons/338b458f-7ebf-449c-9ad1-611eb933b076/concepts/ee21948d-7fad-4821-b61c-0d69bbfcc425

Basically we define the cost function with arguement `fg`.
And the `var` contains all state variables.

And we give `Ipopt` instance the `fg` so that it can store the cost value.

> Student discusses the reasoning behind the chosen N (timestep length) and dt (elapsed duration between timesteps) values. Additionally the student details the previous values tried.

* N = How many iterations we want to predict into future.
* dt =  Length of each iteration in 100 milliseconds.

I started with 25 times into future and dt to be 0.05 like in the lecture. But it was taking too much time to process.
So I made the dt 0.1 but then polyfit sometimes failed due to the fact that a 
3rd polynomial cannot fit into complex curves. So we should keep the future shorter then that. 10 did well enough.
So I predict into 1 second of future in 0.1 precision.

> A polynomial is fitted to waypoints. If the student preprocesses waypoints, the vehicle state, and/or actuators prior to the MPC procedure it is described.

First we register the waypoints in the car's relative perception. So we shift them.
```
for (size_t i = 0; i < ptsx.size(); ++i) {
  const double shift_x = ptsx[i] - px;
  const double shift_y = ptsy[i] - py;
  ptsx[i] = (shift_x * std::cos(0 - psi) -
             shift_y * std::sin(0 - psi));
  ptsy[i] = (shift_x * std::sin(0 - psi) +
             shift_y * std::cos(0 - psi));
}
```

Polyfitting is done in line 112 in main.cpp.
It is 3rd degree polyfitter.
We give it the waypoints, it returns us a vector of 4 coefficients.
```
auto coefficients = polyfit(
    Eigen::Map<Eigen::VectorXd>(ptsx.data(), ptsx.size()),
    Eigen::Map<Eigen::VectorXd>(ptsy.data(), ptsy.size()),
    3);
```

And in line 118 we calculate the initial cross track error compared to origin:
```
double cte = polyeval(coefficients, 0);
```

> The student implements Model Predictive Control that handles a 100 millisecond latency. Student provides details on how they deal with latency.

I tried to compensate a 100ms latency.

For this first I took 2 more state variables from main.cpp line 122,123:

steer_value and throttle value.
```
Eigen::VectorXd state(8);
state << 0, 0, 0, v, cte, epsi, steer_value, throttle_value;
```

It was actually pretty basic, I just increased the starting indices of certain cost predictor loops by one.
Since 1 timestep is 100 ms due to value of dt.

You can see it on mpc.cpp line 202 to 227:
```
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
```

And the main.cpp has the latency embedded in it in line 198:
```
this_thread::sleep_for(chrono::milliseconds(100));
```
