//
// Created by David on 31/01/2022.
//

#ifndef MUJOCO_SANDBOX_ILQR_H
#define MUJOCO_SANDBOX_ILQR_H

#include "../Utility/MujocoController/MujocoController.h"
#include "ilqrCore.h"
#include "iLQR_funcs.h"

using namespace std::chrono;

void iLQR(m_state X0, m_dof *U, m_state *X);
void testILQR();
void loadLastControls(m_dof *U_init);

void  PIDController_Init(PIDController* pid, float Kp, float Ki, float Kd);
float PIDController_Update(PIDController* pid, float setpoint, float measurement);
void myController(const mjModel *m, mjData* d);
void initialseController();

void saveControls(m_dof lastControls, bool fin);
bool poseAchieved(pose desiredPose, pose currentPose);
void initialiseLinearInterpolation(pose _startPose, pose _endPose, float forceMagnitude);

void iLQRSetControlSequence(m_dof *U, int numControls);
void initialiseLinearInterpolation(pose _startPose, pose _endPose, float forceMagnitude);

void setNextControlSequence(const Ref<const VectorXf> U);
void setDesiredRobotConfiguration(const Ref<const m_dof> desiredConfiguration);
void setDesiredEndEffectorPose(pose _desiredEndEffectorPose);

void warmStartControls(m_dof *U, Ref<m_state> X0);

void saveStates(m_state *X_dyn, m_state *X_lin);

#endif //MUJOCO_SANDBOX_ILQR_H
