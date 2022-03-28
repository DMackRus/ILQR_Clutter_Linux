//
// Created by David on 08/02/2022.
//

#ifndef MUJOCO_SANDBOX_ILQR_FUNCS_H
#define MUJOCO_SANDBOX_ILQR_FUNCS_H

#include "../Utility/MujocoController/MujocoController.h"
#include "../Utility/stdInclude/stdInclude.h"


#define MUJOCO_TIMESTEP 0.002
#define DOF         7


void lineariseDynamics(Ref<m_state> currentState, Ref<m_dof> currentControls, Ref<MatrixXd> A, Ref<MatrixXd> B);
void scaleLinearisation(Ref<m_state_state> A_mj_step, Ref<m_state_dof> B_mj_step, Ref<m_state_state> A_dt, Ref<m_state_dof> B_dt, int num_steps_per_dt);
void stepSimulation(const Ref<m_state> currentState, const Ref<m_dof> U, Ref<m_state> Xnew, Ref<m_state> Xdot, int numSimSteps);
float rollOutTrajectory(const Ref<const m_state> X0, m_state *X, m_dof *U, int numControls);
float immediateCost(const Ref<const m_state> X, const Ref<const m_state> X_next, const Ref<const m_dof> U, int controlNum);
float immediateCostAndDerivitives(Ref<m_state> l_x, Ref<m_state_state> l_xx, Ref<m_dof> l_u, Ref<m_dof_dof> l_uu, const Ref<const m_state> X, const Ref<const m_state> X_next, const Ref<const m_dof> U, int controlNum);
m_state costFirstOrderDerivitives(const Ref<const m_state> X, m_state X_next, bool terminal, int controlNum);
float calcStateCost(const Ref<const m_state> X, const Ref<const m_state> X_next, bool terminal, int controlNum);
float calcControlCost(const Ref<const m_dof> U);
float terminalCost(Ref<m_state> l_x, Ref<m_state_state> l_xx, const Ref<const m_state> X);
void initCostMatrices();
void initDesiredState();

void saveTrajecToCSV(m_dof *U, m_state *X);

#endif //MUJOCO_SANDBOX_ILQR_FUNCS_H
