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

#endif //MUJOCO_SANDBOX_ILQR_H
