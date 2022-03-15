//
// Created by davem on 16/02/2022.
//

#include "ilqrCore.h"

extern MujocoController *globalMujocoController;
extern m_state X_desired;
extern float dt;
extern int mujoco_steps_per_dt;

extern int initStateIndex;
extern int baseStateIndex;
extern float horizonLength;
extern int linearising_num_sim_steps;
extern bool alphaSearchEnabled;
extern float maxLamda;
extern float minLamda;
extern float lamdaFactor;
extern float epsConverge;
extern int numControls;
extern int torqueLims[DOF];

extern float angleLims[DOF];
extern float velocityLims[DOF];

float lamb = 0.1;
int numIterations = 0;
float oldCost;

void differentiateDynamics(m_state *X, m_dof *U, m_state_state *f_x, m_state_dof *f_u, float *l, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu){

    MatrixXf I(NUM_STATES, NUM_STATES);
    I.setIdentity(NUM_STATES, NUM_STATES);
    MatrixXf A = ArrayXXf::Zero(NUM_STATES, NUM_STATES);
    MatrixXf B = ArrayXXf::Zero(NUM_STATES, DOF);
    // Linearise the dynamics
    for(int t = 0; t < numControls; t++){
        // Calculate linearised dynamics for current time step via finite differencing
        globalMujocoController->setSystemState(X[t]);
        globalMujocoController->deleteLastMujocoState();
        globalMujocoController->saveMujocoState();
        lineariseDynamics(X[t], U[t], A, B);

//        A *= ((float)linearising_num_sim_steps / mujoco_steps_per_dt);
//        B *= ((float)linearising_num_sim_steps / mujoco_steps_per_dt);

        f_x[t] = A;
        f_u[t] = B;

        l[t] = immediateCostAndDerivitives(l_x[t], l_xx[t], l_u[t], l_uu[t], X[t], X[t+1], U[t]);

        l[t]    *= dt;
        l_x[t]  *= dt;
        l_xx[t] *= dt;
        l_u[t]  *= dt;
        l_uu[t] *= dt;

//        cout << "iteration " << t << endl;
//        cout << " f_x " << f_x[t] << endl;
//        cout << " f_u " << f_u[t] << endl;
//
//        cout << "cost " << l[t] << endl;
//        cout << "cost dif x" << l_x[t] << endl;
//        cout << "cost dif xx" << l_xx[t] << endl;

    }

    l[numControls] = terminalCost(l_x[numControls], l_xx[numControls], X[numControls]);

    l   [numControls] *= dt;
    l_x [numControls] *= dt;
    l_xx[numControls] *= dt;
}

bool backwardsPass(m_state_state *f_x, m_state_dof *f_u, float l, m_state *l_x, m_state_state *l_xx, m_dof *l_u, m_dof_dof *l_uu, m_dof *k,  m_dof_state *K){
    float V = l;
    m_state V_x;
    V_x = l_x[numControls];;
    m_state_state V_xx;
    V_xx = l_xx[numControls];
    m_dof_state f_u_t;
    m_state_state f_x_t;
    bool validBackwardsPass = true;

    for(int t = numControls - 1; t > -1; t--){
        m_state Q_x;
        m_dof Q_u;
        m_state_state Q_xx;
        m_dof_dof Q_uu;
        m_dof_state Q_ux;

        f_u_t = f_u[t].transpose();
        f_x_t = f_x[t].transpose();
        Q_x = l_x[t] + (f_x_t * V_x);
        Q_u = l_u[t] + (f_u_t * V_x);

        Q_xx = l_xx[t] + (f_x_t * (V_xx * f_x[t]));
        Q_ux = (f_u_t * (V_xx * f_x[t]));
        Q_uu = l_uu[t] + (f_u_t * (V_xx * f_u[t]));

        // Caluclate Q_uu_inverse via eigen vector regularisation
        SelfAdjointEigenSolver<m_dof_dof> eigensolver(Q_uu);

        m_dof eigVals = eigensolver.eigenvalues();
        m_dof_dof eigVecs = eigensolver.eigenvectors();

        for(int i = 0; i < DOF; i++){
            if(eigVals(i) < 0) eigVals(i) = 0.0f;
            eigVals(i) += lamb;
            eigVals(i) = 1 / eigVals(i);
        }

        m_dof_dof diagMat;

        diagMat = eigVals.asDiagonal();

        m_dof_dof Q_uu_inv;
        Q_uu_inv = eigVecs * diagMat * eigVecs.transpose();


        k[t] = -Q_uu_inv * Q_u;
        if(isnan(k[t](0))){
            validBackwardsPass = false;
            break;
        }
        K[t] = -Q_uu_inv * Q_ux;

        if(t > 0){
//            cout << "---------------------------------------" << endl;
//            cout << "control number  " << t << endl;
//            cout << "k[t] " << k[t] << endl;
//            cout << "K[t]" << K[t] << endl;
//            cout << "Q_ux" << Q_ux<< endl;
//            cout << "V_xx" << V_xx<< endl;
//            cout << "f_u_t" << f_u_t<< endl;
//            cout << "f_x" << f_x[t]<< endl;
        }

        V_x = Q_x - (K[t].transpose() * (Q_uu * k[t]));
        V_xx = Q_xx - (K[t].transpose() * (Q_uu * K[t]));

    }

    return validBackwardsPass;
}

float forwardsPass(m_state *X, m_state *X_new, m_dof *U, m_dof *U_new, m_dof *k, m_dof_state *K){
    auto forwardPassStart = high_resolution_clock::now();
    bool alphaLineSearch = true;
    float alpha = 1.0f;
    float bestAlphaCost = 0.0f;
    float newAlphaCost = 0.0f;
    bool firstAlpha = true;
    int alphaSearchCount = 0;
    float newCost;
    float initialAlphaCost = 0.0f;

    while(alphaLineSearch){
        m_state xnew(NUM_STATES);
        xnew = X[0];
        globalMujocoController->loadSimulationState(initStateIndex);
        // Calculate new controls using optimal gain matrices
        for(int t = 0; t < numControls; t++){
            m_state stateFeedback;
            stateFeedback = xnew - X[t];
            m_dof feedBackGain = K[t] * stateFeedback;
            m_dof linearFeedback = (alpha * k[t]);

            U_new[t] = U[t] + (alpha * k[t]) + feedBackGain;

            if(t < 5){
//                cout << "---------------------------------------" << endl;
//                cout << "x new " << xnew(7) << endl;
//                cout << "stateFeedback" << stateFeedback(7) << endl;
//                cout << "control number  " << t << endl;
//                cout << "feedbackgain " << feedBackGain(0) << endl;
//                cout << "linearFeedback " << linearFeedback(0) << endl;
//                cout << "old control 1 " << U[t](0) << endl;
//                cout << "new control 1 " << U_new[t](0) << endl;
            }


            // simulate new optimal control and calculate new predicted state.
            m_state _x = xnew;
            m_state _;
            for(int k = 0; k < DOF; k++){
                if(U_new[t](k) > torqueLims[k]) U_new[t](k) = torqueLims[k];
                if(U_new[t](k) < -torqueLims[k]) U_new[t](k) = -torqueLims[k];
            }
            stepSimulation(_x, U_new[t], xnew, _, mujoco_steps_per_dt);
        }

        newAlphaCost = rollOutTrajectory(X[0], X_new, U_new, numControls);

        if(firstAlpha){
            firstAlpha = false;

            bestAlphaCost = newAlphaCost;
            initialAlphaCost = newAlphaCost;
            if(!alphaSearchEnabled){
                break;
            }
            alpha /= 2;
        }
        else{
            float costGrad = abs(newAlphaCost - bestAlphaCost);
            if(costGrad < 0.1){
                alphaLineSearch = false;
            }
            if(newAlphaCost < bestAlphaCost){
                alpha = alpha - pow(0.5, alphaSearchCount + 1);
                bestAlphaCost = newAlphaCost;

            }
            else{
                alpha = alpha + pow(0.5, alphaSearchCount + 1);
            }

        }
        alphaSearchCount++;
    }

    newCost = bestAlphaCost;

    auto forwardPassStop = high_resolution_clock::now();
    auto forwardPassDur = duration_cast<microseconds>(forwardPassStop - forwardPassStart);
    float costImprovement = initialAlphaCost - newCost;

    cout << "Forwards pass: " << forwardPassDur.count()/1000 << " milliseconds" << "num iterations: " << alphaSearchCount << "cost improvement: " << costImprovement << endl;
    return newCost;
}

bool checkForConvergence(float newCost, m_state *X, m_state *X_new, m_dof *U, m_dof *U_new, bool *costImprovement){
    bool convergence = false;
    if(newCost < oldCost){
        *costImprovement = true;
        if(lamb > minLamda){
            lamb /= lamdaFactor;
        }

        // replace old controls with new controls
        // replace old states with new states
        for(int k = 0; k < numControls; k++){
            U[k] = U_new[k];
//            cout << "new control timestep:  " << k << U[k] << endl;
            X[k + 1] = X_new[k + 1];
        }

        std::cout << "--------------------------------------------------" <<  std::endl;
        std::cout << "New cost: " << newCost <<  std::endl;
        std::cout << "Terminal Cube Pos diff " << std::endl;
        m_state X_diff = X[numControls] - X_desired;
        std::cout << "x pos diff: " << X_diff(14) << " z pos diff: " << X_diff(15) << std::endl;

        numIterations++;
        float costGrad = (oldCost - newCost)/newCost;

        if((numIterations > 1) && costGrad < epsConverge){
            convergence = true;
            cout << "ilQR converged, num Iterations: " << numIterations << " final cost: " << newCost << endl;
        }
        oldCost = newCost;
    }
    else{
        // if cost improvement was not made, regulate lamda and run backwards pass again
        lamb *= lamdaFactor;
        if(lamb > maxLamda){
            convergence = true;
            cout << "ilQR finished due to lamda exceeds lamda max " << endl;
        }
    }

    return convergence;
}

void increaseLamda(){
    lamb *= lamdaFactor;
}
