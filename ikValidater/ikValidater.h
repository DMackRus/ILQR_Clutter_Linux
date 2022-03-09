#include "../ikfast/ikfast.h"
#include "../Utility/stdInclude/stdInclude.h"
#include "../Utility/MujocoController/MujocoController.h"
#include "../Utility/MujocoController/MujocoUI.h"

#ifndef IKVALIDATER
#define IKVALIDATER

#define NUM_JOINTS  7
#define PI			3.1415926
#define IKFAST_HAS_LIBRARY

using namespace ikfast;

bool attemptFindValidIKSolution(struct pose endEffector);
bool computeIKSolutions(struct pose endEffector, float freeVal);
m_dof returnValidSolution(int solNumber);
int returnNumValidSols();

bool jointsInsideLimits(m_dof testConfig);


#endif

