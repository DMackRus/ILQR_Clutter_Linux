//
// Created by davem on 21/02/2022.
//

#ifndef MUJOCO_SANDBOX_STDINCLUDE_H
#define MUJOCO_SANDBOX_STDINCLUDE_H

#include <iostream>
#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>

#define NUM_JOINTS 7
#define NUM_STATES 18
// Theta 1 - 7, theta dot 1 - 7, x, z of cube

using namespace Eigen;
using namespace std;
using namespace std::chrono;

struct point3D {
    double x;
    double y;
    double z;
};

struct pose {
    struct point3D pos;
    double roll;
    double pitch;
    double yaw;
};

struct config {
    double jointAngles[NUM_JOINTS];
};

typedef Vector<double, NUM_JOINTS> m_dof;
typedef Vector<double, NUM_STATES> m_state;
typedef Vector<double, 6> m_pose;

typedef Matrix<double, NUM_STATES, NUM_STATES> m_state_state;
typedef Matrix<double, NUM_STATES, NUM_JOINTS> m_state_dof;
typedef Matrix<double, NUM_JOINTS, NUM_STATES> m_dof_state;
typedef Matrix<double, NUM_JOINTS, NUM_JOINTS> m_dof_dof;

#endif //MUJOCO_SANDBOX_STDINCLUDE_H
