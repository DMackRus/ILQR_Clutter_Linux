#include "Utility/MujocoController/MujocoUI.h"
#include "Utility/stdInclude/stdInclude.h"
#include "ikValidater/ikValidater.h"
#include "iLQR/iLQR.h"

#define TRY_IK 0

extern MujocoController *globalMujocoController;

int main() {
    std::cout << "Hello, World!" << std::endl;

    initMujoco();
    globalMujocoController->controlState = 1000;

    pose desiredStartPose;
    desiredStartPose.pos.x = 0.3;
    desiredStartPose.pos.y = 0.3;
    desiredStartPose.pos.z = 0.5;
    desiredStartPose.roll = 0;
    desiredStartPose.pitch = PI;
    desiredStartPose.yaw = 0;

    m_dof initState;
    if(TRY_IK){
        if(attemptFindValidIKSolution(desiredStartPose)){
            initState = returnValidSolution(returnNumValidSols()/2);
        }
        else{
            std::cout << "no valid solution computed for start pose" << std::endl;
            //initState <<
        }
    }
    else{
        // No Y componennt, hovering over table
        initState << -0.564, -0.678, 0.445, -2.45, 0.297, 0.242, -0.297;

        // corner of the table EE facing downwards
        //initState << 0, -0.33, 0.445, -2.2, 0.237, 0.395, 1.93;

        //initState << 0, -0.33, 0.445, -2.2, 0.237, 0.953, 1.93;
    }

    globalMujocoController->setRobotConfiguration(initState);

    testILQR();

//    m_dof desiredState;
//    desiredState << 1.34, 0.825, -1.6, -1.74, 0, 0.165, 0;
//    globalMujocoController->setDesiredRobotConfiguration(desiredState);
    //globalMujocoController->controlState = linearInterpolation;

//    pose desiredEE;
//    desiredEE.pos.x = 0.7;
//    desiredEE.pos.y = 0;
//    desiredEE.pos.z = 0.45;
//    desiredEE.pos.x = 0.4;
//    desiredEE.pos.y = 0;
//    desiredEE.pos.z = 0.8;
    //globalMujocoController->setDesiredEndEffectorPose(desiredEE);

    //pose startPose = globalMujocoController->returnEndEffectorPos();

    //initialiseLinearInterpolation(startPose, desiredEE, 20);


    render();



    return 0;
}
