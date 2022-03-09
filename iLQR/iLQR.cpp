//
// Created by David on 31/01/2022.
//

#include "iLQR.h"
// find CHECK, FIX, FILL

#define LOAD_LAST_CONTROLS 1

extern MujocoController *globalMujocoController;

extern int numControls;
extern float oldCost;
extern int mujoco_steps_per_dt;

void iLQR(m_state X0, m_dof *U, m_state *X){
    bool optimisationFinished = false;
    int numIterations = 0;
    float newCost = 0;

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
    globalMujocoController->controlState = ilqrSim;
    VectorXf X0(NUM_STATES);
    X0 << -0.564, -0.678, 0.445, -2.45, 0.297, 0.242, -0.297,
            0, 0, 0, 0, 0, 0, 0,
            0.41, 0;

    initCostMatrices();
    initDesiredState();

    m_dof *initControls = new m_dof[numControls];
    if(loadLastControls){
        loadLastControls(initControls);
    }
    else{
        warmStartControls(initControls, X0);
    }

    for(int i = 0; i < numControls; i++){
        cout << "control joint 0, index " << i << " " << initControls[i](0) << endl;
    }


    m_state *X_lin = new m_state[numControls + 1];
    m_state *X_dyn = new m_state[numControls + 1];
    X_lin[0] = X0;
    X_dyn[0] = X0;

    auto setSimStateStart = high_resolution_clock::now();
    globalMujocoController->setSystemState(X0);

    auto setSimStateStop = high_resolution_clock::now();
    auto setStateDur = duration_cast<microseconds>(setSimStateStop - setSimStateStart);

    cout << "setting simulation state " << setStateDur.count() << " microseconds" << endl;
    globalMujocoController->saveMujocoState();
    globalMujocoController->saveMujocoState();
    globalMujocoController->loadSimulationState(0);

//    MatrixXf A = ArrayXXf::Zero(14, 14);
//    MatrixXf B = ArrayXXf::Zero(14, 7);
//    m_state _;
//    for(int i = 0;  i < numControls; i++){
//        globalMujocoController->setSystemState(X_dyn[i]);
//        globalMujocoController->deleteLastMujocoState();
//        globalMujocoController->saveMujocoState();
//        lineariseDynamics(X_dyn[i], initControls[i], A, B);
//        stepSimulation(X_dyn[i], initControls[i], X_dyn[i+1], _, mujoco_steps_per_dt);
//
////        cout << "------------------------- A ------------------- " << endl;
////        cout << A << endl;
////        cout << "------------------------- B ------------------- " << endl;
////        cout << B << endl;
//        A *= ((float)linearising_num_sim_steps / mujoco_steps_per_dt);
//        B *= ((float)linearising_num_sim_steps / mujoco_steps_per_dt);
//        cout << "------------------------- A ------------------- " << endl;
//        cout << A << endl;
//        cout << "------------------------- B ------------------- " << endl;
//        cout << B << endl;
//
//        m_state xdot_lin;
//        m_state actual_x_dot;
//        xdot_lin = (((A * X_lin[i]) + (B * initControls[i])) * dt);
//        actual_x_dot = X_dyn[i+1] - X_dyn[i];
//        cout << "xdot was: " << endl;
//        cout << xdot_lin << endl;
//        cout << "actual x dot was: " << endl;
//        cout << actual_x_dot << endl;
//        X_lin[i+1] = X_lin[i] + xdot_lin;
//        //        cout << "------------------------- B ------------------- " << endl;
////        cout << B << endl;
//
//
//    }
//
//    m_state diff;
//    float sumDiff[7] = {0, 0, 0, 0, 0, 0, 0};
//    for(int i = 0; i < 7; i++){
//        diff(i) = X_dyn[numControls - 1](i) - X_lin[numControls - 1](i);
//        sumDiff[i] += pow(diff(i), 2);
//    }
//    cout << "sum squared diff at end: " << endl;
//    for(int i = 0; i < 7; i++){
//        cout << sumDiff[i] << " "  << endl;
//    }
//    saveStates(X_dyn, X_lin);

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
    globalMujocoController->iLQRSetControlSequence(initControls, numControls);
    globalMujocoController->controlState = simulating;
    globalMujocoController->mujocoTimesStepsPerControl = mujoco_steps_per_dt;
    globalMujocoController->numControlsPerTrajectory = numControls;
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