#include "Utility/MujocoController/MujocoUI.h"
#include "iLQR/iLQR.h"
#include "Utility/stdInclude/stdInclude.h"

#define ILQR 1

extern MujocoController *globalMujocoController;
extern int controlState;



// Eigen is row, col
int main() {

    initMujoco();
    initialseController();
    // Starting joint angles - 2.25, 0.825, -1.66, -1.81, 0, -0.258, 0
    // cube starting near corner of the table - position -
//    double epsilon = 0.1;
//    double startpos = 1;
//    double startVel = 0.1;
//    double force = 1;
//    double posStates[4] = {0};
//    double velStates[4] = {0};
//    posStates[0] = startpos;
//    velStates[0] = startVel;
//    double incPosStates[4] = {0};
//    double incVelStates[4] = {0};
//    incPosStates[0] = startpos;
//    incVelStates[0] = startVel + epsilon;
//    double dt = 0.1;
//
//    for(int i = 0; i < 3; i++){
//        double accel = force - (0.2 * velStates[i]);
//        velStates[i+1] = velStates[i] + (accel * dt);
//        posStates[i+1] = posStates[i] + (velStates[i+1] * dt);
//
//        accel = force - (0.2 * incVelStates[i]);
//        incVelStates[i+1] = incVelStates[i] + (accel * dt);
//        incPosStates[i+1] = incPosStates[i] + (incVelStates[i+1] * dt);
//    }
//
//    double posDiff = (incPosStates[3] - posStates[3]) / epsilon;
//    double velDiff = (incVelStates[3] - velStates[3]) / epsilon;
//    int a = 1;
    if(ILQR){
        testILQR();
    }
    else{
        simpleTest();
    }




    render();

    return 1;

}

