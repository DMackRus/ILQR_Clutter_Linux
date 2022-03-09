#include "ikValidater.h"
#include "../ikfast/ikfast.h"

m_dof lastSolFound;
m_dof *validSolutions = new m_dof[200];
bool printOnce = true;
extern MujocoController* globalMujocoController;
int solutionIndex = 0;

using namespace ikfast;

bool attemptFindValidIKSolution(struct pose endEffector) {
	solutionIndex = 0;
	float lowerVFreeLim = -2.97;
	float upperVFreeLim = 2.97;
	int numSteps = 100;
	float diff = (upperVFreeLim - lowerVFreeLim) / (numSteps);
	int counter = 0;
	bool validSolutionFound = false;
	bool computing = true;
	float vFree = lowerVFreeLim;
	while (computing == true) {
		// check a solution for values of vfree

		if (computeIKSolutions(endEffector, vFree)) {
			validSolutionFound = true;
			validSolutions[solutionIndex] = lastSolFound;
			solutionIndex++;
		}
		vFree += diff;

		// if valid solution found exit out of while loop
		counter++;
		// can exit out of the loop if a valid solution is found, or compute all of them
		if (validSolutionFound == true) {
			//computing = false;
		}
		if (counter > numSteps) {
			computing = false;
		}
	}

	return validSolutionFound;
}

bool computeIKSolutions(struct pose endEffector, float freeVal) {
	IkSolutionList<IkReal> solutions;
	std::vector<IkReal> vfree(GetNumFreeParameters());
	IkReal eeRot[9], eeTrans[3];

	//printf("roll: %f, pitch: %f, yaw: %f \n", pose[0], pose[1], pose[2]);

	double su = sin(endEffector.roll);
	double cu = cos(endEffector.roll);
	double sv = sin(endEffector.pitch);
	double cv = cos(endEffector.pitch);
	double sw = sin(endEffector.yaw);
	double cw = cos(endEffector.yaw);

	// top row of rot matrix
	eeRot[0] = cv * cw;
	eeRot[1] = (su * sv * cw) - (cu * sw);
	eeRot[2] = (su * sw) + (cu * sv * cw);

	// second row of rot matrix
	eeRot[3] = cv * sw;
	eeRot[4] = (cu * cw) + (su * sv * sw);
	eeRot[5] = (cu * sv * sw) - (su * cw);

	// third row of rot matrix
	eeRot[6] = -sv;
	eeRot[7] = su * cv;
	eeRot[8] = cu * cv;

	if (printOnce) {
		printOnce = false;
		printf("roll: %f, pitch: %f, yaw: %f \n", endEffector.roll, endEffector.pitch, endEffector.yaw);
		printf("[ %f %f %f \n ", eeRot[0], eeRot[1], eeRot[2]);
		printf("  %f %f %f \n ", eeRot[3], eeRot[4], eeRot[5]);
		printf("  %f %f %f ] \n ", eeRot[6], eeRot[7], eeRot[8]);
	}


	eeTrans[0] = endEffector.pos.x;
	eeTrans[1] = endEffector.pos.y;
	eeTrans[2] = endEffector.pos.z;

	for (std::size_t i = 0; i < vfree.size(); i++) {
		vfree[i] = freeVal;
	}

	bool bSuccess = ComputeIk(eeTrans, eeRot, &vfree[0], solutions);

	if (!bSuccess) {
		return false;
	}

	//printf("Found %d ik solutions:\n", (int)solutions.GetNumSolutions());
	std::vector<IkReal> solvalues(GetNumJoints());
	
	// Only the last two solutions seem valid from the computeIk Function
	for (int solutionIndex = solutions.GetNumSolutions() - 2; solutionIndex < solutions.GetNumSolutions(); solutionIndex++) {
		const IkSolutionBase<IkReal>& sol = solutions.GetSolution(solutionIndex);
		std::vector<IkReal> vsolfree(sol.GetFree().size());
		sol.GetSolution(&solvalues[0], vsolfree.size() > 0 ? &vsolfree[0] : NULL);
		float ikSolution[NUM_JOINTS];

		for (char i = 0; i < NUM_JOINTS; i++) {
			if (i == 5) {
				if (solvalues[i] > 0) {
					ikSolution[i] = solvalues[i] - (PI / 2);
				}
				else {
					ikSolution[i] = solvalues[i] + PI + (PI / 2);
				}
			}
			else {
				ikSolution[i] = solvalues[i];
			}
			lastSolFound(i) = ikSolution[i];
		}
		
		// if greater than 0 subtract pi / 4
		if (lastSolFound(6) > 0) {
			lastSolFound(6) = lastSolFound(6) - (PI / 4);
		}
		// else add 3 pi by 4
		else {
			lastSolFound(6) = lastSolFound(6) + (3 * PI / 4);
		}


		

		globalMujocoController->setRobotConfiguration(lastSolFound);
        updateScreen();

		if (!(globalMujocoController->isConfigInCollision(lastSolFound))) {
            if (jointsInsideLimits(lastSolFound)) {
                return true;
            }
            //return true;
		}
        else{
            int a = 1;
        }

	}

	/*printf("(free=%d): ", (int)sol.GetFree().size());
	for (int i = 0; i < NUM_JOINTS; ++i)
		printf("%.15f, ", ikSolution[i]);
	printf("\n");*/

	return false;
}

m_dof returnValidSolution(int solNumber) {
	return validSolutions[solNumber];
}

int returnNumValidSols() {
	return solutionIndex;
}

// Inverse kinematics sometimes specifies joint angles thatr are outside the operating range of the joint
bool jointsInsideLimits(m_dof testConfig) {
	bool jointsValid = true;

	if (testConfig(0) > 2.97 || testConfig(0) < -2.97) {
		jointsValid = false;
        std::cout << "joint 0 failed: " << testConfig(0) << std::endl;
	}

	// Hard joint limit
	if (testConfig(1) > 1.83 || testConfig(1) < -1.83) {
		jointsValid = false;
        std::cout << "joint 1 failed: " << testConfig(1) << std::endl;
	}

	if (testConfig(2) > 2.97 || testConfig(2) < -2.97) {
		jointsValid = false;
        std::cout << "joint 2 failed: " << testConfig(2) << std::endl;
	}

	if (testConfig(3) > 3.14 || testConfig(3) < 0) {
		jointsValid = false;
        std::cout << "joint 3 failed: " << testConfig(3) << std::endl;
	}

	if (testConfig(4) > 2.97 || testConfig(4) < -2.97) {
		jointsValid = false;
        std::cout << "joint 4 failed: " << testConfig(4) << std::endl;
	}

	// Hard joint limit 
	if (testConfig(5) > 2.18 || testConfig(5) < -1.66) {
		jointsValid = false;
        std::cout << "joint 5 failed: " << testConfig(5) << std::endl;
	}

	return jointsValid;
}


