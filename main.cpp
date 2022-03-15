#include "Utility/MujocoController/MujocoUI.h"
#include "iLQR/iLQR.h"
#include "Utility/stdInclude/stdInclude.h"

#define ILQR 1

extern MujocoController *globalMujocoController;



// Eigen is row, col
int main() {

    initMujoco();
    initialseController();
    // Starting joint angles - 2.25, 0.825, -1.66, -1.81, 0, -0.258, 0
    // cube starting near corner of the table - position -

//    m_dof initState;
//    initState << -0.564, -0.678, 0.445, -2.45, 0.297, 0.242, -0.297;
//    globalMujocoController->setRobotConfiguration(initState);


    if(ILQR){
        testILQR();
    }

    render();

    return 1;

}

