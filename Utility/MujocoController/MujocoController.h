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
    m_dof* controlSequence;
    m_dof nextControlSequence;
    int controlCounter = 0;
    int numControlsPerTrajectory;
    int mujocoTimesStepsPerControl;
    int mujocoTimeStepCounter = 0;
    bool resetSimFlag = false;
    int controlState = 0;
    PIDController positionPID[NUM_JOINTS];
    m_dof desiredJointAngles;

    pose desiredEndEffectorPos;
    pose startEndEffectorPos;
    pose diffPoseStartDesired;
    Vector3d linearInterpolationDesiredForce;

    MatrixXd J_COMi;// = MatrixXd::Zero(3, _model->nv); //COM Jacobian of body i with shape 3x number of degrees of freedom

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

    void setDesiredRobotConfiguration(const Ref<const m_dof> desiredConfiguration);
    void setDesiredEndEffectorPose(pose _desiredEndEffectorPose);
    struct pose returnEndEffectorPos();

    bool isConfigInCollision(m_dof configuration);
    int getRobotNumberOfCollisions();

    void saveSimulationState();
    void setResetLevel(int resetLevel);
    void resetSimulation();
    void loadSimulationState(int stateIndex);
    void resetSimulationToStart();


    void setNextControlSequence(const Ref<const VectorXf> U);
    void step();
    void iLQRSetControlSequence(m_dof *U, int numControls);
    void initialiseLinearInterpolation(pose _startPose, pose _endPose, float forceMagnitude);


    Eigen::MatrixXd calculateJacobian(int bodyId);
};

void  PIDController_Init(PIDController* pid, float Kp, float Ki, float Kd);
float PIDController_Update(PIDController* pid, float setpoint, float measurement);
void myController(const mjModel *m, mjData* d);
void initialseController();

void saveControls(m_dof lastControls, bool fin);
bool poseAchieved(pose desiredPose, pose currentPose);
void initialiseLinearInterpolation(pose _startPose, pose _endPose, float forceMagnitude);

#endif //CLIONS_PROJECTS_MUJOCOCONTROLLER_H
