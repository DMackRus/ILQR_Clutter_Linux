//
// Created by David on 31/01/2022.
//

#include "iLQR.h"
// find CHECK, FIX, FILL
//extern mjfGeneric mjcb_control;
int controlState = 0;
PIDController positionPID[NUM_JOINTS];
m_dof desiredJointAngles;

m_dof* controlSequence;
m_dof lastControl;
int controlCounter = 0;
int numControlsPerTrajectory;
int mujocoTimesStepsPerControl;
int mujocoTimeStepCounter = 0;

pose desiredEndEffectorPos;
pose startEndEffectorPos;
pose diffPoseStartDesired;
Vector3d linearInterpolationDesiredForce;

extern MujocoController *globalMujocoController;

extern int numControls;
extern float oldCost;
extern int mujoco_steps_per_dt;

extern float torqueLims[7];
extern m_state X_desired;
extern m_dof nextControlSequence;

extern int linearising_num_sim_steps;
extern float dt;

ofstream outputDiffDyn;

std::string diffDynFilename = "diffDyn.csv";

#define TEST_LINEARISATION 1

void iLQR(m_state X0, m_dof *U, m_state *X){
    bool optimisationFinished = false;
    int numIterations = 0;
    float newCost = 0;
    controlState = ilqrSim;

    // Initialise partial differentiation matrices for all timesteps T
    // for linearised dynamics
    m_state_state *f_x = new m_state_state[numControls];
    m_state_dof *f_u = new m_state_dof[numControls];

    // Quadratic cost partial derivitives
    float l[numControls];
    m_state *l_x = new m_state[numControls + 1];
    m_state_state *l_xx = new m_state_state[numControls + 1];
    m_dof *l_u = new m_dof[numControls];
    m_dof_dof *l_uu = new m_dof_dof[numControls];

    // Initialise state feedback gain matrices
    m_dof *k = new m_dof[numControls];
    m_dof_state *K = new m_dof_state[numControls];

    // Initialise new controls and states storage for evaluation
    m_dof *U_new = new m_dof[numControls];
    m_state *X_new = new m_state[numControls + 1];

    globalMujocoController->saveMujocoState();

    //Perform a forward pass with initialised controls -FILL
    oldCost = rollOutTrajectory(X0, X, U, numControls);

    // iterate until optimisation finished
    while(!optimisationFinished){

        // Linearise the dynamics and save cost values at each state
        auto start = high_resolution_clock::now();

        // STEP 1 - Linearise dynamics and calculate cost quadratics at every time step
        differentiateDynamics(X, U, f_x, f_u, l, l_x, l_xx, l_u, l_uu);

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Linearising model: " << duration.count()/1000 << " milliseconds" << endl;

        bool costImprovementMade = false;
        while(!costImprovementMade){

            // STEP 2 - Backwards pass to compute optimal linear and feedback gain matrices k and K
            if(!backwardsPass(f_x, f_u, l[numControls - 1], l_x, l_xx, l_u, l_uu, k, K)){
                increaseLamda();
            }
            else{
                // STEP 3 - Forwards pass to calculate new optimal controls - with optional alpha binary search
                newCost = forwardsPass(X, X_new, U, U_new, k, K);

                // STEP 4 - Check for convergence
                optimisationFinished = checkForConvergence(newCost, X, X_new, U, U_new, &costImprovementMade);

                if(optimisationFinished == true){
                    break;
                }
            }
        }
    }
}

void testILQR(){
    controlState = ilqrSim;
    VectorXf X0(NUM_STATES);
    X0 << -0.564, -0.678, 0.445, -2.45, 0.297, 0.242, -0.297,
            0, 0, 0, 0, 0, 0, 0,
            0.41, 0;
    //0.41 start for cube

    initCostMatrices();
    initDesiredState();

    globalMujocoController->setSystemState(X0);

    globalMujocoController->saveMujocoState();
    globalMujocoController->saveMujocoState();
    globalMujocoController->loadSimulationState(0);

    m_dof *initControls = new m_dof[numControls];

    warmStartControls(initControls, X0);

//    for(int i = numControls/2; i < numControls; i++){
//        for(int j = 0; j < DOF; j++){
//            initControls[i](j) = 0;
//        }
//    }

    m_state *X_lin = new m_state[numControls + 1];
    m_state *X_dyn = new m_state[numControls + 1];
    X_lin[0] = X0;
    X_dyn[0] = X0;



    globalMujocoController->loadSimulationState(0);
    controlState = ilqrSim;

    MatrixXf A = ArrayXXf::Zero(NUM_STATES, NUM_STATES);
    MatrixXf B = ArrayXXf::Zero(NUM_STATES, DOF);
    m_state _;

//    for(int i = 0; i < 5; i ++){
//        globalMujocoController->setSystemState(X_dyn[i]);
//        stepSimulation(X_dyn[i], initControls[0], X_dyn[i+1], _, 1);
//
//
//    }

    m_state diff;
    float sumDiff[16] = {0};

    if(TEST_LINEARISATION){
        for(int i = 0;  i < numControls; i++){
            //globalMujocoController->setSystemState(X_dyn[i]);
            globalMujocoController->deleteLastMujocoState();
            globalMujocoController->saveMujocoState();
            lineariseDynamics(X_dyn[i], initControls[i], A, B);
            stepSimulation(X_dyn[i], initControls[i], X_dyn[i+1], _, mujoco_steps_per_dt);

//        cout << "------------------- x dynamic ---------------------" << endl;
//        cout << X_dyn[i+1] << endl;
//        A *= ((float)linearising_num_sim_steps / mujoco_steps_per_dt);
//        B *= ((float)linearising_num_sim_steps / mujoco_steps_per_dt);
//        cout << "------------------------- A ------------------- " << endl;
//        cout << A << endl;
//        cout << "------------------------- B ------------------- " << endl;
//        cout << B << endl;

            m_state xdot_lin;
            m_state actual_x_dot;
            X_lin[i+1] = ((A * X_dyn[i]) + (B * initControls[i]));
            actual_x_dot = X_dyn[i+1] - X_dyn[i];
            xdot_lin = X_lin[i+1] - X_lin[i];
//        cout << "xdot was: " << endl;
//        cout << xdot_lin << endl;
//        cout << "actual x dot was: " << endl;
//        cout << actual_x_dot << endl;
            //X_lin[i+1] = X_lin[i] + xdot_lin;
            diff = X_dyn[i] - X_lin[i];
            for(int  j = 0; j < 16; j++){
                sumDiff[j] += pow((X_dyn[i](j) - X_lin[i](j)), 2);

            }
            int a = 1;

        }

        cout << "sum squared diff at end: " << endl;
        float totalBadness = 0.0f;
        for(int i = 0; i < 16; i++){
            cout << sumDiff[i] << " "  << endl;
            totalBadness += sumDiff[i];
        }
        cout << "total badness: " << totalBadness << endl;
        saveStates(X_dyn, X_lin);
    }

    globalMujocoController->loadSimulationState(0);
    auto iLQRStart = high_resolution_clock::now();

    iLQR(X0, initControls, X_dyn);

    auto iLQRStop = high_resolution_clock::now();
    auto iLQRDur = duration_cast<microseconds>(iLQRStop - iLQRStart);

    cout << "iLQR took " << iLQRDur.count()/1000000 << " seconds" << endl;

    saveTrajecToCSV(initControls, X_dyn);

    // reset simulation
    // save control sequence to mujoco controller class
    // untick ilqractive
    globalMujocoController->loadSimulationState(0);
    iLQRSetControlSequence(initControls, numControls);
    controlState = simulating;
    mujocoTimesStepsPerControl = mujoco_steps_per_dt;
    numControlsPerTrajectory = numControls;
}

void loadLastControls(m_dof *U_init){
    ifstream in("initControls.csv");
    vector<vector<double>> fields;

    if (in) {
        string line;

        int linecounter = 0;
        while (getline(in, line)) {
            stringstream sep(line);
            string field;

            int counter = 0;
            while (getline(sep, field, ',')) {
                U_init[linecounter](counter) = stof(field);
                counter++;
            }
            linecounter++;
        }
    }
    int a = 1;
}

void  PIDController_Init(PIDController* pid, float Kp, float Ki, float Kd) {
    pid->integrator = 0.0f;
    pid->prevError = 0.0f;
    pid->differentiator = 0.0f;
    pid->prevMeasurement = 0.0f;
    pid->out = 0.0f;

    pid->Kp = Kp;
    pid->Ki = Ki;
    pid->Kd = Kd;

    pid->tau = 0.002f;
    pid->T = 0.002f;

    pid->limMinInt = -100.0f;
    pid->limMaxInt = 100.0f;

    pid->limMax = 87.0f;
    pid->limMin = -87.0f;
}

float PIDController_Update(PIDController* pid, float setpoint, float measurement) {
    float error = setpoint - measurement;


    /*
    * Proportional
    */
    float proportional = pid->Kp * error;


    /*
    * Integral
    */
    pid->integrator = pid->integrator + (pid->Ki * pid->T * ((error + pid->prevError) / 2));

    /* Anti-wind-up via integrator clamping */
    if (pid->integrator > pid->limMaxInt) {

        pid->integrator = pid->limMaxInt;

    }
    else if (pid->integrator < pid->limMinInt) {

        pid->integrator = pid->limMinInt;

    }

    /*
    * Derivative (band-limited differentiator)
    */

    pid->differentiator = pid->Kd * ((measurement - pid->prevMeasurement) / pid->T);

    //pid->differentiator = -(2.0f * pid->Kd * (measurement - pid->prevMeasurement)	/* Note: derivative on measurement, therefore minus sign in front of equation! */
    //	+ (2.0f * pid->tau - pid->T) * pid->differentiator)
    //	/ (2.0f * pid->tau + pid->T);


    /*
    * Compute output and apply limits
    */
    pid->out = proportional + pid->integrator + pid->differentiator;
    //pid->out = proportional + pid->integrator;
    //pid->out = proportional;

    if (pid->out > pid->limMax) {

        pid->out = pid->limMax;

    }
    else if (pid->out < pid->limMin) {

        pid->out = pid->limMin;

    }

    /* Store error and measurement for later use */
    pid->prevError = error;
    pid->prevMeasurement = measurement;

    /* Return controller output */
    return pid->out;
}

//void saveControls(m_dof lastControls, bool fin){
//    saveControlsFile << lastControls(0) << "," << lastControls(1) << "," << lastControls(2) << "," << lastControls(3) << "," << lastControls(4) << "," << lastControls(5) << "," << lastControls(6) << endl;
//    if(fin){
//        saveControlsFile.close();
//    }
//}

void initialiseLinearInterpolation(pose _startPose, pose _endPose, float forceMagnitude){
    startEndEffectorPos = _startPose;
    pose diff;

    diff.pos.x = desiredEndEffectorPos.pos.x - startEndEffectorPos.pos.x;
    diff.pos.y = desiredEndEffectorPos.pos.y - startEndEffectorPos.pos.y;
    diff.pos.z = desiredEndEffectorPos.pos.z - startEndEffectorPos.pos.z;
    diffPoseStartDesired = diff;
    float magnitudeDiff = sqrt(pow(diff.pos.x,2) + pow(diff.pos.y,2) + pow(diff.pos.z,2));

    diff.pos.x /= magnitudeDiff;
    diff.pos.y /= magnitudeDiff;
    diff.pos.z /= magnitudeDiff;

    linearInterpolationDesiredForce(0) = diff.pos.x * forceMagnitude;
    linearInterpolationDesiredForce(1) = diff.pos.y * forceMagnitude;
    linearInterpolationDesiredForce(2) = diff.pos.z * forceMagnitude;
    cout << "linear interpolation desired force: x:" << linearInterpolationDesiredForce(0) << " y: " << linearInterpolationDesiredForce(1) << " z: " << linearInterpolationDesiredForce(2) << endl;
}

bool poseAchieved(pose desiredPose, pose currentPose){
    bool poseAcheived = false;
    float diff = sqrt(pow((desiredPose.pos.x - currentPose.pos.x),2)
                      + pow((desiredPose.pos.y - currentPose.pos.y),2)
                      + pow((desiredPose.pos.z - currentPose.pos.z),2));

    if(diff < 0.1){
        poseAcheived = true;
    }
    return poseAcheived;
}

void myController(const mjModel *m, mjData *d){
    m_dof nextControl;

    if(controlState == ilqrSim){
        nextControl = nextControlSequence;
        for(int i = 0; i < NUM_JOINTS; i++) {
            d->ctrl[i] = nextControlSequence(i);
        }
    }

    else if(controlState == simulating){
        nextControl = controlSequence[controlCounter];

        for(int i = 0; i < NUM_JOINTS; i++){
            d->ctrl[i] = controlSequence[controlCounter](i);

        }

        mujocoTimeStepCounter++;
        if(mujocoTimeStepCounter >= mujocoTimesStepsPerControl){
            m_state currentState = globalMujocoController->returnSystemState();
            mujocoTimeStepCounter = 0;
            controlCounter++;
        }

        if(controlCounter >= numControlsPerTrajectory){
            controlCounter = 0;
            globalMujocoController->resetSimFlag = true;
            m_state endState = globalMujocoController->returnSystemState();
        }
    }

    else if(controlState == staticPos){
        double measurments[NUM_JOINTS];
        float newControls[NUM_JOINTS];

        for (int i = 0; i < NUM_JOINTS; i++) {
            measurments[i] = d->qpos[i];
        }

        for (int i = 0; i < NUM_JOINTS; i++) {
            newControls[i] = PIDController_Update(&positionPID[i], desiredJointAngles(i), measurments[i]);
            nextControl(i) = newControls[i];
        }

        for (int i = 0; i < NUM_JOINTS; i++) {
            d->ctrl[i] = newControls[i];
        }


    }
    // Linearly interpolate end effector between two 3D positions whilst countering gravity
    else if(controlState == linearInterpolation){
        m_dof testSaveControls;
        int endEffectorId = globalMujocoController->returnModelID("franka_gripper");

        MatrixXd jacobian = globalMujocoController->calculateJacobian(endEffectorId);
        MatrixXd jacobian_T = jacobian.transpose();
        pose currentEndEffector = globalMujocoController->returnEndEffectorPos();

        float endEffecPosDiff[3];
        endEffecPosDiff[0] = currentEndEffector.pos.x - startEndEffectorPos.pos.x;
        endEffecPosDiff[1] = currentEndEffector.pos.y - startEndEffectorPos.pos.y;
        endEffecPosDiff[2] = currentEndEffector.pos.z - startEndEffectorPos.pos.z;

        float percentageAchieved[3];
        percentageAchieved[0] = (endEffecPosDiff[0] / diffPoseStartDesired.pos.x) * 100;
        percentageAchieved[1] = (endEffecPosDiff[1] / diffPoseStartDesired.pos.y) * 100;
        percentageAchieved[2] = (endEffecPosDiff[2] / diffPoseStartDesired.pos.z) * 100;

        Vector3d correctiveEndEffecForce;

        correctiveEndEffecForce(0) = (linearInterpolationDesiredForce(0));
        correctiveEndEffecForce(1) = (linearInterpolationDesiredForce(1));
        correctiveEndEffecForce(2) = (linearInterpolationDesiredForce(2));

        // if the end effector is below correct height
        float redFactor = 100 * endEffecPosDiff[2];
        correctiveEndEffecForce(2) -= redFactor;

        VectorXd desiredJointTorques = jacobian_T * correctiveEndEffecForce;

        for( int i = 0; i < NUM_JOINTS; i++){
            d->ctrl[i] = d->qfrc_bias[i] + desiredJointTorques(i);
            nextControl(i) = d->ctrl[i];
        }

        if(poseAchieved(desiredEndEffectorPos, currentEndEffector)){
            controlState = staticCalc;
        }
    }
    else if(controlState == staticCalc){
        for( int i = 0; i < NUM_JOINTS; i++){
            d->ctrl[i] = d->qfrc_bias[i];
            nextControl(i) = d->ctrl[i];
        }
    }
    else{

    }
    lastControl = nextControl;
}

void iLQRSetControlSequence(m_dof *U, int numControls){

    controlSequence = new m_dof[numControls];
    for(int i = 0; i < numControls; i++){
        controlSequence[i] = U[i];
    }
}

void initialseController(){
    mjcb_control = myController;
    PIDController_Init(&positionPID[0], 870, 10, 0.5);
    PIDController_Init(&positionPID[1], 870, 10, 0.5);
    PIDController_Init(&positionPID[2], 870, 10, 0.5);
    PIDController_Init(&positionPID[3], 870, 10, 0.5);
    PIDController_Init(&positionPID[4], 120, 10, 0.5);
    PIDController_Init(&positionPID[5], 120, 10, 0.5);
    PIDController_Init(&positionPID[6], 120, 10, 0.5);
}

void warmStartControls(m_dof *U, Ref<m_state> X0){
    // set the initial configuration of the robot arm

    // use a Jacobian linear interpolation controller to move forwards at a set pace
    controlState = linearInterpolation;

    pose desiredEE;
    desiredEE.pos.x = 0.7;
    desiredEE.pos.y = 0;
    desiredEE.pos.z = 0.45;
    setDesiredEndEffectorPose(desiredEE);
    pose startPose = globalMujocoController->returnEndEffectorPos();

    initialiseLinearInterpolation(startPose, desiredEE, 40);

    for(int i = 0; i < numControls; i++){
        for(int j = 0; j < mujoco_steps_per_dt; j++){
            globalMujocoController->step();
        }
        U[i] = lastControl;
    }
}

void setDesiredRobotConfiguration(const Ref<const m_dof> desiredConfiguration){
    desiredJointAngles = desiredConfiguration;
}

void setDesiredEndEffectorPose(pose _desiredEndEffectorPose){
    desiredEndEffectorPos = _desiredEndEffectorPose;
}

void saveStates(m_state *X_dyn, m_state *X_lin){
    ofstream outputDiffDyn;

    std::string diffDynFilename = "diffDyn.csv";
    outputDiffDyn.open(diffDynFilename);
    outputDiffDyn << "X0 dyn" << "," << "X0 lin" << "," << "X0 diff" << "," << "X1 dyn" << "," << "X1 lin" << "," << "X1 diff" << ",";
    outputDiffDyn << "X2 dyn" << "," << "X2 lin" << "," << "X2 diff" << "," << "X3 dyn" << "," << "X3 lin" << "," << "X3 diff" << ",";
    outputDiffDyn << "X4 dyn" << "," << "X4 lin" << "," << "X4 diff" << "," << "X5 dyn" << "," << "X5 lin" << "," << "X5 diff" << ",";
    outputDiffDyn << "X6 dyn" << "," << "X6 lin" << "," << "X6 diff" << "," << "X7 dyn" << "," << "X7 lin" << "," << "X7 diff" << ",";
    outputDiffDyn << "X8 dyn" << "," << "X8 lin" << "," << "X8 diff" << "," << "X9 dyn" << "," << "X9 lin" << "," << "X9 diff" << ",";
    outputDiffDyn << "X10 dyn" << "," << "X10 lin" << "," << "X10 diff" << "," << "X11 dyn" << "," << "X11 lin" << "," << "X11 diff" << ",";
    outputDiffDyn << "X12 dyn" << "," << "X12 lin" << "," << "X12 diff" << "," << "X13 dyn" << "," << "X13 lin" << "," << "X13 diff" << ",";
    outputDiffDyn << "X14 dyn" << "," << "X14 lin" << "," << "X14 diff" << "," << "X15 dyn" << "," << "X15 lin" << "," << "X15 diff" << "," << endl;

    for(int i = 0; i < numControls; i++){
        for(int j = 0; j < NUM_STATES; j++){
            outputDiffDyn << X_dyn[i](j) << ",";
            outputDiffDyn << X_lin[i](j) << ",";
            outputDiffDyn << X_lin[i](j) - X_dyn[i](j) << ",";
        }
        outputDiffDyn << endl;
    }

}