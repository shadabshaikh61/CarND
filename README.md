# CarND-Controls-MPC
Self-Driving Car Engineer Nanodegree Program

---

## Reflections

### Intro

This contains my solution for MPC project of Udacity SDC nano degree program. Our goal is to navigate car through whole track in simulator without getting off the road. This is controlled mainly by changing accelaration and steering control.

The solution for this problem uses C++ libraries namely Ipopt and CppAD for calculating optimum path of vehicle. Also the path is modeled as 3rd-degree polynomial. Cost function for vehicle is calculated by taaking into consideration velocity, cross track error, orientation, etc of the vehicle.

### Implementation

#### Model

This is the actual kinematics of the motion model. This motion model take into account vehicle x & y coordinates, car's orientation(psi), steering angle and distance from center of road which is cross-track-error(cte). By considering all the mentioned parameters kinematics models and two controls parameters namely accelaration and steering control next state of car can be predicted with help of previous state from equation shown below.

![equations](./eqns.png)

#### Timestep Length and Elapsed Duration (N & dt)

The value i have choosed for the is project for N and dt are 10 & 0.1, respectively. N is number of points and dt is time differences between two points that need to be consider for cost calculation. If dt is choosen very large, it is very likely that model will predict the output which is highly incorrect. And this may lead to overshoot of vehicle while low value of dt lead to prediction of very small distance. If N is choosen very large then more computation is required for cost calculation and if N is choosen very small the very less points are calculated which will not be able to give correct output. I have tried different N and dt such that prediction of position of next 1 second can be found.

#### Polynomial Fitting and MPC Preprocessing

The waypoints are transformed such that those point 

``x_new = x_old * cos(psi) + y_old * sin(psi)``

``y_new = y_old * cos(psi) - x_old * sin(psi)``

#### Model Predictive Control with Latency - 100ms

---

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1(mac, linux), 3.81(Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets
    cd uWebSockets
    git checkout e94b6e1
    ```
    Some function signatures have changed in v0.14.x. See [this PR](https://github.com/udacity/CarND-MPC-Project/pull/3) for more details.

* **Ipopt and CppAD:** Please refer to [this document](https://github.com/udacity/CarND-MPC-Project/blob/master/install_Ipopt_CppAD.md) for installation instructions.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). This is already part of the repo so you shouldn't have to worry about it.
* Simulator. You can download these from the [releases tab](https://github.com/udacity/self-driving-car-sim/releases).
* Not a dependency but read the [DATA.md](./DATA.md) for a description of the data sent back from the simulator.


## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./mpc`.
