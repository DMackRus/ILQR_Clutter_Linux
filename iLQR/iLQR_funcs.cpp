//
// Created by David on 08/02/2022.
//

#include "iLQR_funcs.h"

/**************************************************************************
 *
 *  iLQR Parameters
 *
 *
 */
float horizonLength = 15.0; // seconds
float dt = 0.04; // time between controls changing
int numControls = horizonLength / dt; // number of controls needed as per horizon length and time step
int mujoco_steps_per_dt = (dt / MUJOCO_TIMESTEP) + 1; // Number of mj_steps needed per control time step
int linearising_num_sim_steps = 5;  // How many mujoco steps to use for linearising
bool alphaSearchEnabled = false;     // Whether alpha search is enabled to maximise optimisation at each forwards pass
float maxLamda = 10000;             // Maximum lambda before canceliing optimisation
float minLamda = 0.00001;            // Minimum lamda
float lamdaFactor = 10;             // Lamda multiplicative factor
float epsConverge = 0.005;
bool oneSidedFiniteDiff = true;
bool costFunctionFD = false;
m_dof nextControlSequence;


//float controlCost[DOF] = {0.0001, 0.0001, 0.0001, 0.0001, 0.00005, 0.00005, 0.00005};
float controlCost[DOF] = {0, 0, 0, 0, 0, 0, 0};
float stateCosts[NUM_STATES] = {0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0,
                                    0.5, 2, 0.1, 0.1};
float terminalScalarConstant = 5;

/**************************************************************************
 *
 *  Other variables
 *
 *
 */

int torqueLims[DOF] = {87, 87, 87, 87, 12, 12, 12};
float angleLims[DOF] = {2.9, 1.76, 2.89, 1.5, 2.9, 1.89, 2.9};
float velocityLims[DOF] = {2.175, 2.175, 2.175, 2.175, 2.6, 2.6, 2.6};

int baseStateIndex = 1;
int initStateIndex = 0;

m_state X_desired;
m_dof_dof R;
m_state_state Q;
m_state_state Q_term;
m_dof_dof Z;

ofstream outputFile;

std::string filename = "iLQR.csv";


extern MujocoController *globalMujocoController;


// Make a linear model of the dynamics of the system at the current state
// Create an approximation of the type x(t.) = Ax(t) + Bu(t)
// This function calculates the A and B matrices using finite differencing
void lineariseDynamics(Ref<m_state> currentState, Ref<m_dof> currentControls, Ref<MatrixXd> A, Ref<MatrixXd> B){

//    float epsState[NUM_STATES] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1};
    float epsState[NUM_STATES] = {1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6};
    float epsControls[DOF] = {1e-3,1e-3,1e-3,0.5,0.5,1e-3,0.5};
    float eps_state = 1e-3;
    float eps_control = 1e-3;
    float eps = 1e-6;
    A.setIdentity();

    auto linDynStart = high_resolution_clock::now();
    int microSecs_Loading = 0;
    auto startTimerLoading = high_resolution_clock::now();
    auto stopTimerLoading = high_resolution_clock::now();
    auto loadingDuration = duration_cast<microseconds>(stopTimerLoading - startTimerLoading);
    //int numberMujocoStepsNeeded = linearising_num_sim_steps;
    int numberMujocoStepsNeeded = 1;

    // Make a copy of the current state
    m_state X_copy = currentState.replicate(1, 1);

    // calculate the A matrix
    for(int i = DOF - 1; i < NUM_STATES; i++){
        // Create an incremented and decremented version of current state
        if((i != 14) && (i != 15)){
            m_state X_inc = X_copy.replicate(1,1);
            X_inc(i) += eps_state;
            m_state X_dec = X_copy.replicate(1, 1);
            X_dec(i) -= eps_state;

            // apply same controls and observe how state variables change with respect to changing individual state variables
            m_state stateInc(NUM_STATES);
            m_state stateDec(NUM_STATES);
            m_state _inc(NUM_STATES);
            m_state _dec(NUM_STATES);
            m_state normalStateAfterStep;
//        startTimerLoading = high_resolution_clock::now();
//        globalMujocoController->loadSimulationState(baseStateIndex);
//        stepSimulation(X_copy, currentControls, normalStateAfterStep, stateInc, numberMujocoStepsNeeded);
//        cout << "----------------- state after step no increment --------------------- " << endl;
//        cout << normalStateAfterStep << endl;
            globalMujocoController->loadSimulationState(baseStateIndex);
            globalMujocoController->setSystemState(X_inc);

//        stopTimerLoading = high_resolution_clock::now();
//        loadingDuration = duration_cast<microseconds>(stopTimerLoading - startTimerLoading);
//        microSecs_Loading += loadingDuration.count();
            stepSimulation(X_inc, currentControls, _inc, stateInc, numberMujocoStepsNeeded);
            //stepSimulation(X_copy, currentControls, _inc, stateInc, numberMujocoStepsNeeded);
//        cout << "----------------- state increment --------------------- " << endl;
//        cout << _inc << endl;
//        startTimerLoading = high_resolution_clock::now();

            globalMujocoController->loadSimulationState(baseStateIndex);
            globalMujocoController->setSystemState(X_dec);

//        stopTimerLoading = high_resolution_clock::now();
//        loadingDuration = duration_cast<microseconds>(stopTimerLoading - startTimerLoading);
//        microSecs_Loading += loadingDuration.count();
            stepSimulation(X_dec, currentControls, _dec, stateDec, numberMujocoStepsNeeded);
//        cout << "----------------- state decrement --------------------- " << endl;
//        cout << _dec << endl;

            // calculate the gradient
            for(int j = 0; j < NUM_STATES; j++){
                // CHECK, needs double checking especially with mujoco
                A(j, i) = (_inc(j) - _dec(j) ) / (2 * eps_state);
            }

//        cout << "------------------- A so far ---------------------" << endl;
//        cout << A << endl;
            int a = 1;
        }
    }

//    for(int i = 0; i < NUM_STATES; i++) {
//        for (int j = 0; j < NUM_STATES; j++) {
//            if(A(i, j) > 5){
//
//                cout << A << endl;
//                int a = 1;
//            }
//        }
//    }
//
    //reduce veloicty rows dependant on position to zero for stability
    for(int i = 0; i < NUM_STATES; i++){
        for(int j = 0; j < NUM_STATES; j++){
            if(i <= 7){
                if(j >= 14 && j <= 15){
                    A(i, j) = 0;
                }
            }


            if(i >= 7 && i <= 13){
                if(j <= 6){
                    A(i, j) = 0;
                }
                if(j == 14 || j == 15) {
                    A(i, j) = 0;
                }
            }
            if(i >= 14 && i <=15){
                if(j <= 6){
                    A(i, j) = 0;
                }
            }

            if(i >= 16){
                if(j <= 6){
                    A(i, j) = 0;
                }
                if(j == 14 || j == 15) {
                    A(i, j) = 0;
                }
            }


        }
    }

//    cout << "------------------- A so far ---------------------" << endl;
//    cout << A << endl;


    m_dof U_copy = currentControls.replicate(1, 1);
    for(int i = 0; i < DOF; i++){
        m_dof U_inc = U_copy.replicate(1, 1);
        m_dof U_dec = U_copy.replicate(1, 1);

        U_inc(i) += eps_control;
        U_dec(i) -= eps_control;

        m_state stateInc(NUM_STATES);
        m_state stateDec(NUM_STATES);
        m_state _inc(NUM_STATES);
        m_state _dec(NUM_STATES);

//        startTimerLoading = high_resolution_clock::now();

        globalMujocoController->loadSimulationState(baseStateIndex);

//        stopTimerLoading = high_resolution_clock::now();
//        loadingDuration = duration_cast<microseconds>(stopTimerLoading - startTimerLoading);
//        microSecs_Loading += loadingDuration.count();
        stepSimulation(X_copy, U_inc, _inc, stateInc, numberMujocoStepsNeeded);
//        cout << "----------------- x dot increment --------------------- " << endl;
//        cout << _inc << endl;

//        startTimerLoading = high_resolution_clock::now();

        globalMujocoController->loadSimulationState(baseStateIndex);
        //globalMujocoController->setSystemState(X_copy, accels);

//        stopTimerLoading = high_resolution_clock::now();
//        loadingDuration = duration_cast<microseconds>(stopTimerLoading - startTimerLoading);
//        microSecs_Loading += loadingDuration.count();

        stepSimulation(X_copy, U_dec, _dec, stateDec, numberMujocoStepsNeeded);
//        cout << "----------------- x dot decrement --------------------- " << endl;
//        cout << _dec << endl;

        for(int j = 0; j < NUM_STATES; j++){
            B(j, i) = (_inc(j) - _dec(j))/(2 * eps_control);
        }

    }
//    cout << "------------------- B so far ---------------------" << endl;
//    cout << B << endl;

    globalMujocoController->loadSimulationState(baseStateIndex);

//    auto linDynStop = high_resolution_clock::now();
//    auto linDynDur = duration_cast<microseconds>(linDynStop - linDynStart);

//    cout << "Time taken by linearising dynamics: " << linDynDur.count() << " microseconds" << endl;
//    cout << "Time taken loading was " << microSecs_Loading << endl;
}

void scaleLinearisation(Ref<m_state_state> A_mj_step, Ref<m_state_dof> B_mj_step, Ref<m_state_state> A_dt, Ref<m_state_dof> B_dt, int num_steps_per_dt){
    m_state_dof lastBTerm = B_mj_step;
    B_dt = lastBTerm;
    A_dt = A_mj_step;
    for(int j = 0; j < num_steps_per_dt - 1; j++){
        lastBTerm = A_mj_step * lastBTerm;
        B_dt += lastBTerm;
        A_dt = A_dt * A_mj_step;
    }
}

void stepSimulation(const Ref<m_state> currentState, const Ref<m_dof> U, Ref<m_state> Xnew, Ref<m_state> Xdot, int numSimSteps){

    m_dof acc = globalMujocoController->returnRobotAccelerations();
    m_dof vel = globalMujocoController->returnRobotVelocities();
    m_dof pos = globalMujocoController->returnRobotConfiguration();

//    cout << "robot acc before:" << endl;
//    cout << acc << endl;
//    cout << "robot vel before:" << endl;
//    cout << vel << endl;
//    cout << "robot pos before:" << endl;
//    cout << pos << endl;

    nextControlSequence = U;
    for(int i = 0; i < numSimSteps; i++){
        globalMujocoController->step();
    }

    Xnew = globalMujocoController->returnSystemState();
    acc = globalMujocoController->returnRobotAccelerations();
    vel = globalMujocoController->returnRobotVelocities();
    pos = globalMujocoController->returnRobotConfiguration();

//    cout << "old state" << endl;
//    cout << currentState << endl;
//    cout << "new state " << endl;
//    cout << Xnew << endl;

//    cout << "robot acc after:" << endl;
//    cout << acc << endl;
//    cout << "robot vel after:" << endl;
//    cout << vel << endl;
//    cout << "robot pos after:" << endl;
//    cout << pos << endl;

    m_state x_dot_differently;

    for(int i = 0; i < NUM_STATES; i++) {
        Xdot(i) = (Xnew(i) - currentState(i)) / (MUJOCO_TIMESTEP * numSimSteps);

    }

    for(int i = 0; i < DOF; i++){
        Xdot(i) = vel(i);
        Xdot(i + 7) = acc(i);
    }
//    cout << "x dot is" << endl;
//    cout << Xdot << endl;
//    cout << "x dot alternate method" << endl;
//    cout << x_dot_differently << endl;
    int a = 1;

}

float rollOutTrajectory(const Ref<const m_state> X0, m_state *X, m_dof *U, int numControls){

    float cost = 0.0f;
    float l;
    globalMujocoController->loadSimulationState(initStateIndex);
    X[0] = X0;
    m_state l_x;
    m_state_state l_xx;
    m_dof l_u;
    m_dof_dof l_uu;

    for(int i = 0; i < numControls; i++){
        // Calculate cost associated with current state Xt and current control Ut

//        cout << "Current control in rollout, iteration: " << i << " " << U[i] << endl;
//        cout << "Current state in rollout, iteration: " << i << " " << X[i] << endl;
        l = immediateCost(X[i], X[i+1], U[i], i);
        cost += (l * dt);

        // Step simulation set number of times
        nextControlSequence = U[i];
        for(int j = 0; j < mujoco_steps_per_dt; j++){
            globalMujocoController->step();
        }


        // update new state variable Xt+1
        X[i+1] = globalMujocoController->returnSystemState();

    }
    return cost;
}

float immediateCost(const Ref<const m_state> X, const Ref<const m_state> X_next, const Ref<const m_dof> U, int controlNum){
    float cost;
    float eps = 1e-1;
    m_state X_diff;

    // actual - desired
    X_diff = X - X_desired;

    cost = calcStateCost(X, X_next, false, controlNum) + calcControlCost(U);

    return cost;
}

float immediateCostAndDerivitives(Ref<m_state> l_x, Ref<m_state_state> l_xx, Ref<m_dof> l_u, Ref<m_dof_dof> l_uu, const Ref<const m_state> X, const Ref<const m_state> X_next, const Ref<const m_dof> U, int controlNum){
    float cost;
    float eps = 1e-1;
    m_state X_diff;

    // actual - desired
    X_diff = X - X_desired;

    cost = calcStateCost(X, X_next, false, controlNum) + calcControlCost(U);

    if(costFunctionFD){
        l_x = costFirstOrderDerivitives(X, X_next, false, controlNum);
        for(int i = 0; i < NUM_STATES; i++){
            m_state X_inc = X.replicate(1, 1);
            m_state X_dec = X.replicate(1, 1);

            X_inc(i) += eps;
            X_dec(i) -= eps;

            m_state l_x_inc = costFirstOrderDerivitives(X_inc, X_next, false, controlNum);
//        cout << "l_x_inc" << endl;
//        cout << l_x_inc << endl;
            m_state l_x_dec = costFirstOrderDerivitives(X_dec, X_next, false, controlNum);
//        cout << "l_x_dec" << endl;
//        cout << l_x_dec << endl;

            for(int j = 0 ; j < NUM_STATES; j++){
                l_xx(j, i) = (l_x_inc(j) - l_x_dec(j)) / (2 * eps);
            }
        }
    }
    else{
        float power = (float) controlNum / numControls;
        float scalar = pow(terminalScalarConstant, power);
        l_x = 2 * scalar * Q * X_diff;
        l_xx = 2 * scalar * Q;
    }

//    cout << "X_diff was" << endl;
//    cout << X_diff << endl;
//
//    m_state orig_l_x = (Q * X_diff);
//    m_state_state orig_l_xx = Q;
//    cout << "original l_x" << endl;
//    cout << orig_l_x << endl;
//
//    cout << "new l_x" << endl;
//    cout << l_x << endl;
//
//    cout << "original l_xx" << endl;
//    cout << orig_l_xx << endl;
//
//    cout << "new l_xx" << endl;
//    cout << l_xx << endl;

    l_u = R * U;
    l_uu = R;

    return cost;
}

float terminalCost(Ref<m_state> l_x, Ref<m_state_state> l_xx, const Ref<const m_state> X){
    float cost;
    float eps = 1e-1;

    m_state X_diff = X - X_desired;
    cost = calcStateCost(X, X, true, numControls);

    if(costFunctionFD){
        l_x = costFirstOrderDerivitives(X, X, true, 0);

        for(int i = 0; i < NUM_STATES; i++){
            m_state X_inc = X.replicate(1, 1);
            m_state X_dec = X.replicate(1, 1);

            X_inc(i) += eps;
            X_dec(i) -= eps;

            m_state l_x_inc = costFirstOrderDerivitives(X_inc, X_inc, true, 0);
            m_state l_x_dec = costFirstOrderDerivitives(X_dec, X_inc, true, 0);

            for(int j = 0 ; j < NUM_STATES; j++){
                l_xx(j, i) = (l_x_inc(j) - l_x_dec(j)) / (2 * eps);
            }
        }
    }
    else{
        l_x = X_diff.transpose() * Q_term;
        l_xx = Q_term;
    }

    return cost;
}


m_state costFirstOrderDerivitives(const Ref<const m_state> X, m_state X_next, bool terminal, int controlNum){
    m_state l_x;
    m_state X_inc;
    m_state X_dec;
    float eps = 1e-2;

    for(int i = 0; i < NUM_STATES; i++){
        X_inc = X.replicate(1, 1);
        X_dec = X.replicate(1, 1);

        X_inc(i) += eps;
        X_dec(i) -= eps;

        float incCost = calcStateCost(X_inc, X_next, terminal, controlNum);
        float decCost = calcStateCost(X_dec, X_next, terminal, controlNum);

        l_x(i) = (incCost - decCost) / (2 * eps);
    }


    return l_x;
}

float calcStateCost(const Ref<const m_state> X,const Ref<const m_state> X_next, bool terminal, int controlNum){
    float stateCost;
    m_state X_diff;
    m_dof accel;
    VectorXd temp(1);
    VectorXd tempAccel(1);

    // actual - desired
    X_diff = X - X_desired;
    float power = controlNum / numControls;
    float scalar = pow(terminalScalarConstant, power);

    if(terminal){
        temp = scalar * X_diff.transpose() * Q_term * X_diff;
    }
    else{
        temp = scalar * X_diff.transpose() * Q * X_diff;
    }

    stateCost = temp(0);


    return stateCost;
}

float calcControlCost(const Ref<const m_dof> U){
    float controlCost;
    VectorXd temp(1);

    temp = 0.5 * U.transpose() * R * U;
    controlCost = temp(0);

    return controlCost;
}



void initCostMatrices(){

    R.setIdentity();
    for(int i = 0; i < DOF; i++){
        R(i, i) = controlCost[i];
    }
    Q.setIdentity();
    for(int i = 0; i < NUM_STATES; i++){
        Q(i, i) = stateCosts[i];
    }
    Q_term = Q.replicate(1, 1);
    for(int i = 0; i < NUM_STATES; i++){
        Q_term(i, i) *= 10;
    }
}

void initDesiredState(){
    X_desired << 0, 0, 0, -0.4, 0, 0, 0,
                 0, 0, 0, 0, 0, 0, 0,
                 0.7, -0.2, 0, 0;

    // Currently its the last two entries in the desired state that matter, which is the x and y position of the cube

    outputFile.open(filename);
    outputFile << "Control Number" << "," << "T1" << "," << "T2" << "," << "T3" << "," << "T4" << "," << "T5" << "," << "T6" << "," << "T7" << ",";
    outputFile << "V1" << "," << "V2" << "," << "V3" << "," << "V4" << "," << "V5" << "," << "V6" << "," << "V7" << ",";
    outputFile << "p1" << "," << "p2" << "," << "p3" << "," << "p4" << "," << "p5" << "," << "p6" << "," << "p7" << endl;

}

void saveTrajecToCSV(m_dof *U, m_state *X){
    for(int i = 0; i < numControls; i++){
        outputFile << i << "," << U[i](0) << "," << U[i](1) << "," << U[i](2) << "," << U[i](3) << "," << U[i](4) << "," << U[i](5) << "," << U[i](6) << ",";

        outputFile << X[i](7);
        for(int j = 0; j < DOF - 1; j++){
            outputFile << "," << X[i](j+1+7);
        }

        outputFile << "," << X[i](0);
        for(int j = 0; j < DOF - 1; j++){
            outputFile << "," << X[i](j+1);
        }
        outputFile << endl;
    }
    outputFile.close();
}