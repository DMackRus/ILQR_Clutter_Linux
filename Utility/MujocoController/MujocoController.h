//
// Created by davem on 20/01/2022.
//
#ifndef CLIONS_PROJECTS_MUJOCOCONTROLLER_H
#define CLIONS_PROJECTS_MUJOCOCONTROLLER_H

#include "mujoco.h"
#include "glfw3.h"
#include "../stdInclude/stdInclude.h"


using namespace std;

#define NUM_JOINTS  7

typedef struct {

    /* Controller gains */
    float Kp;
    float Ki;
    float Kd;

    /* Derivative low-pass filter time constant */
    float tau;

    /* Output limits */
    float limMin;
    float limMax;

    /* Integrator limits */
    float limMinInt;
    float limMaxInt;

    /* Sample time (in seconds) */
    float T;

    /* Controller "memory" */
    float integrator;
    float prevError;			/* Required for integrator */
    float differentiator;
    float prevMeasurement;		/* Required for differentiator */

    /* Controller output */
    float out;

} PIDController;

enum controlStates{
    simulating,
    ilqrSim,
    staticPos,
    linearInterpolation,
    staticCalc
};

class MujocoController {
public:
    MujocoController(mjModel* m, mjData* d);

    mjModel* _model;
    mjData* _data;
    int robotBodyID[11];
    int _resetLevel = 0;
    std::vector<mjData*> _mujocoStates;

    bool resetSimFlag = false;

    MatrixXd J_COMi;

    mjModel* returnModel();
    mjData* returnCurrentModelData();

    // Save current mujoco state and add it to a array of mujoco states
    void saveMujocoState();

    // delete last entry in array of mujoco states
    void deleteLastMujocoState();

    // Directly altering or returning the state of the simulation (position, velocity, acceleration)
    void setSystemState(const Ref<const m_state> systemState);
    m_state returnSystemState();
    void setBodyState(int bodyId, const Ref<const m_pose> pose);
    m_pose returnBodyState(int bodyId);
    void setRobotConfiguration(const Ref<const VectorXf> configuration);
    ArrayXf returnRobotConfiguration();
    void setRobotVelocities(const Ref<const VectorXf> jointVelocities);
    ArrayXf returnRobotVelocities();
    void setRobotAccelerations(const Ref<const VectorXf> jointAccelerations);
    ArrayXf returnRobotAccelerations();

    struct pose returnEndEffectorPos();

    bool isConfigInCollision(m_dof configuration);
    int getRobotNumberOfCollisions();

    void saveSimulationState();
    void setResetLevel(int resetLevel);
    void resetSimulation();
    void loadSimulationState(int stateIndex);
    void resetSimulationToStart();

    void step();

    Eigen::MatrixXd calculateJacobian(int bodyId);
    int returnModelID(const std::string& input);
};



#endif //CLIONS_PROJECTS_MUJOCOCONTROLLER_H