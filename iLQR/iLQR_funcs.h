//
// Created by David on 08/02/2022.
//

#ifndef MUJOCO_SANDBOX_ILQR_FUNCS_H
#define MUJOCO_SANDBOX_ILQR_FUNCS_H

#include "../Utility/MujocoController/MujocoController.h"
#include "../Utility/stdInclude/stdInclude.h"


#define MUJOCO_TIMESTEP 0.002
#define DOF         7


void lineariseDynamics(Ref<m_state> currentState, Ref<m_dof> currentControls, Ref<MatrixXf> A, Ref<MatrixXf> B);
void stepSimulation(const Ref<m_state> currentState, const Ref<m_dof> U, Ref<m_state> Xnew, Ref<m_state> Xdot, int numSimSteps);
float rollOutTrajectory(const Ref<const VectorXf> X0, m_state *X, m_dof *U, int numControls);
float immediateCost(const Ref<const m_state> X, const Ref<const m_state> X_next, const Ref<const m_dof> U);
float immediateCostAndDerivitives(Ref<VectorXf> l_x, Ref<MatrixXf> l_xx, Ref<VectorXf> l_u, Ref<MatrixXf> l_uu, const Ref<const VectorXf> X, const Ref<const VectorXf> X_next, const Ref<const VectorXf> U);
m_state costFirstOrderDerivitives(const Ref<const m_state> X, m_state X_next, bool terminal);
float calcStateCost(const Ref<const m_state> X, const Ref<const m_state> X_next, bool terminal);
float calcControlCost(const Ref<const m_dof> U);
float terminalCost(Ref<VectorXf> l_x, Ref<MatrixXf> l_xx, const Ref<const VectorXf> X);
void initCostMatrices();
void initDesiredState();

void saveTrajecToCSV(m_dof *U, m_state *X);

#endif //MUJOCO_SANDBOX_ILQR_FUNCS_H
