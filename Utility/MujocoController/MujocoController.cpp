//
// Created by davem on 20/01/2022.
//
#include "MujocoController.h"

mjModel* model;						// MuJoCo model
mjData* mdata;						// MuJoCo data
mjvCamera cam;						// abstract camera
mjvOption opt;						// visualization options
mjvScene scn;						// abstract scene
mjrContext con;						// custom GPU context
mjvPerturb pert;
GLFWwindow* window;
MujocoController* globalMujocoController;

int printCounter = 0;
int controlInitCounter = 0;
m_state testingState;

ofstream saveControlsFile;

std::string saveControlsFilename = "lastControls.csv";
int saveControlCounter = 0;
int numControlsToSaveCounter = 0;

int _numControls = 300;
int _mujoco_steps_per_dt = 10;


MujocoController::MujocoController(mjModel* m, mjData* d){
    _model = m;
    _data = d;

    PIDController_Init(&positionPID[0], 870, 10, 0.5);
    PIDController_Init(&positionPID[1], 870, 10, 0.5);
    PIDController_Init(&positionPID[2], 870, 10, 0.5);
    PIDController_Init(&positionPID[3], 870, 10, 0.5);
    PIDController_Init(&positionPID[4], 120, 10, 0.5);
    PIDController_Init(&positionPID[5], 120, 10, 0.5);
    PIDController_Init(&positionPID[6], 120, 10, 0.5);

    saveControlsFile.open(saveControlsFilename);
}

void MujocoController::step(){
    mj_step(_model, _data);
}

mjModel* MujocoController::returnModel() {
    return _model;
}

//mjData* MujocoController::returnCurrentModelData(){
//
//    return newData;
//}

void MujocoController::saveMujocoState(){
    mjData *newData = mj_makeData(_model);
    mj_copyData(newData, _model, _data);
    _mujocoStates.push_back(newData);

}
void MujocoController::deleteLastMujocoState(){
    mj_deleteData(_mujocoStates[1]);
    _mujocoStates.pop_back();
}

void MujocoController::setSystemState(const Ref<const m_state> systemState){
    m_dof robotConfig;
    m_dof robotVelocities;

    for(int i = 0; i < NUM_JOINTS; i++){
        robotConfig(i) = systemState(i);
        robotVelocities(i) = systemState(i + 7);
        //robotAccelerations(i) = systemState(i + 14);
    }
    setRobotConfiguration(robotConfig);
    setRobotVelocities(robotVelocities);

    int boxId = mj_name2id(model, mjOBJ_BODY, "box_obstacle_1");
    m_pose boxPose = returnBodyState(boxId);

    // Setting the x and y position to the desired state
    boxPose(0) = systemState(14);
    boxPose(1) = systemState(15);

    setBodyState(boxId, boxPose);

}

m_state MujocoController::returnSystemState(){
    m_state systemState;
    for(int i = 0; i < NUM_JOINTS; i++){
        systemState(i) = _data->qpos[i];
        systemState(i + 7) = _data->qvel[i];
        //systemState(i + 14) = _data->qacc[i];
    }

    int boxId = mj_name2id(model, mjOBJ_BODY, "box_obstacle_1");
    systemState(14) = _data->xpos[3 * boxId];
    systemState(15) = _data->xpos[(3 * boxId) + 1];

    return systemState;
}

void MujocoController::setBodyState(int bodyId, const Ref<const m_pose> pose){

    // TODO extend this to work with rotations also, whether its quaternions or euler angles
//    for(int i = 0; i < 3; i++){
//        _data->qpos[(bodyId * 3) + i] = pose(i);
//    }

    _data->qpos[16] = pose(0);
    _data->qpos[17] = pose(1);
}

m_pose MujocoController::returnBodyState(int bodyId){
    m_pose bodyPose;

    // TODO extend this to work with rotations also, whether its quaternions or euler angles
    for(int i = 0; i < 3; i++){
        bodyPose(i) = _data->xpos[(bodyId * 3) + i];
    }

    return bodyPose;
}



void MujocoController::setRobotConfiguration(const Ref<const VectorXf> configuration) {

    for (int i = 0; i < NUM_JOINTS; i++) {
        _data->qpos[i] = configuration(i);
    }
    mj_forward(model, _data);
}

ArrayXf MujocoController::returnRobotConfiguration(){
    ArrayXf robotConfig(7);

    for(int i = 0; i < NUM_JOINTS; i++){
        robotConfig(i) = _data->qpos[i];
    }
    return robotConfig;
}

void MujocoController::setRobotVelocities(const Ref<const VectorXf> jointVelocities){
    for (int i = 0; i < NUM_JOINTS; i++) {
        _data->qvel[i] = jointVelocities(i);
    }
}

ArrayXf MujocoController::returnRobotVelocities(){
    ArrayXf robotVelocities(7);
    for(int i = 0; i < NUM_JOINTS; i++){
        robotVelocities(i) = _data->qvel[i];
    }
    return robotVelocities;
}

void MujocoController::setRobotAccelerations(const Ref<const VectorXf> jointAccelerations){
    for (int i = 0; i < NUM_JOINTS; i++) {
        _data->qacc[i] = jointAccelerations(i);
    }
}

ArrayXf MujocoController::returnRobotAccelerations(){
    ArrayXf jointAccelerations(NUM_JOINTS);
    for(int i = 0; i < NUM_JOINTS; i++){
        jointAccelerations(i) = _data->qacc[i];
    }


    return jointAccelerations;
}

void MujocoController::setDesiredRobotConfiguration(const Ref<const m_dof> desiredConfiguration){
    desiredJointAngles = desiredConfiguration;
}

void MujocoController::setDesiredEndEffectorPose(pose _desiredEndEffectorPose){
    desiredEndEffectorPos = _desiredEndEffectorPose;
}

bool MujocoController::isConfigInCollision(m_dof configuration) {
    bool collision = false;
    int originalResetValue = 0;

    setRobotConfiguration(configuration);
    mj_step1(_model, _data);

    int numberOfCollisions = getRobotNumberOfCollisions();

    mj_resetData(_model, _data);
    if (_resetLevel > 0) {
        mjData* targetData = _mujocoStates.at(_resetLevel - 1);
        mj_copyData(_data, _model, targetData);
    }

    //Appears not to need this
    //mj_forward(m, d);

    if (numberOfCollisions > 0) {
        collision = true;
    }

    return collision;
}

int MujocoController::getRobotNumberOfCollisions() {

    int numContacts = _data->ncon;
    int numCollisions = 0;
    for (int i = 0; i < numContacts; i++) {
        auto contact = _data->contact[i];

        // Get the ids of the two bodies in contacts
        int bodyInContact1 = _model->body_rootid[_model->geom_bodyid[contact.geom1]];
        int bodyInContact2 = _model->body_rootid[_model->geom_bodyid[contact.geom2]];

        // only consider it a collision if robot - robot
        // or robot - table

        bool contact1Robot = false;
        bool contact1Table = false;
        bool contact2Robot = false;
        bool contact2Table = false;
        for (int j = 0; j < 11; j++) {
            if (bodyInContact1 == robotBodyID[j]) {
                contact1Robot = true;
            }

            if (bodyInContact2 == robotBodyID[j]) {
                contact2Robot = true;
            }
        }

//        if (bodyInContact1 == tableId) {
//            contact1Table = true;
//        }
//        if (bodyInContact2 == tableId) {
//            contact2Table = true;
//        }

        if (contact1Robot) {
            if (contact2Robot || contact2Table) {
                numCollisions++;
            }
        }
        else if(contact2Robot) {
            if (contact1Robot || contact1Table) {
                numCollisions++;
            }
        }
    }

    return numCollisions;
}

struct pose MujocoController::returnEndEffectorPos(){
    pose endEffectorPose;

    // FIX add code to assign rotation to endeffector pose
    endEffectorPose.pos.x = _data->xpos[30];
    endEffectorPose.pos.y = _data->xpos[31];
    endEffectorPose.pos.z = _data->xpos[32];

    return endEffectorPose;
}

void MujocoController::saveSimulationState(){
    mjData *newData = mj_makeData(_model);
    mj_copyData(newData, _model, _data);
    _mujocoStates.push_back(newData);

}

void MujocoController::setResetLevel(int resetLevel){
    _resetLevel = resetLevel;
}

void MujocoController::resetSimulation(){
    mj_resetData(_model, _data);

    if (_resetLevel > 0) {
        mjData *targetData = _mujocoStates.at(_resetLevel);
        mj_copyData(_data, _model, targetData);
    }

    mj_forward(_model, _data);
}

void MujocoController::loadSimulationState(int stateIndex){

    if(0){
        mj_copyData(_data, _model, _mujocoStates[stateIndex]);
    }
    else{
        _data->time = _mujocoStates[stateIndex]->time;
        mju_copy(_data->qpos, _mujocoStates[stateIndex]->qpos, _model->nq);
        mju_copy(_data->qvel, _mujocoStates[stateIndex]->qvel, _model->nv);
        mju_copy(_data->act,  _mujocoStates[stateIndex]->act,  _model->na);

        // copy mocap body pose and userdata
        mju_copy(_data->mocap_pos,  _mujocoStates[stateIndex]->mocap_pos,  3*_model->nmocap);
        mju_copy(_data->mocap_quat, _mujocoStates[stateIndex]->mocap_quat, 4*_model->nmocap);
        mju_copy(_data->userdata, _mujocoStates[stateIndex]->userdata, _model->nuserdata);

        // copy warm-start acceleration
        mju_copy(_data->qacc_warmstart, _mujocoStates[stateIndex]->qacc_warmstart, _model->nv);
    }
}

Eigen::MatrixXd MujocoController::calculateJacobian(int bodyId){
    Eigen::MatrixXd kinematicJacobian(3, 7);

    mjtNum* J_COMi_temp = mj_stackAlloc(_data, 3*_model->nv);

    mj_jacBody(_model, _data, J_COMi_temp, NULL, bodyId);

    J_COMi = Map<Matrix<double, Dynamic, Dynamic, RowMajor>>(J_COMi_temp, 3, _model->nv);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 7; j++) {
            kinematicJacobian(i, j) = J_COMi(i, j);
            //cout << kinematicJacobian(i, j) << endl;
        }
    }

    return kinematicJacobian;
}

void MujocoController::setNextControlSequence(const Ref<const VectorXf> U){
    nextControlSequence = U;
}

void myController(const mjModel *m, mjData *d){


    if(globalMujocoController->controlState == ilqrSim){
        for(int i = 0; i < NUM_JOINTS; i++) {
            d->ctrl[i] = globalMujocoController->nextControlSequence(i);
        }
    }

    else if(globalMujocoController->controlState == simulating){
        for(int i = 0; i < NUM_JOINTS; i++){
            d->ctrl[i] = globalMujocoController->controlSequence[globalMujocoController->controlCounter](i);
        }
        globalMujocoController->mujocoTimeStepCounter++;
        if(globalMujocoController->mujocoTimeStepCounter >= globalMujocoController->mujocoTimesStepsPerControl){
            globalMujocoController->mujocoTimeStepCounter = 0;
            globalMujocoController->controlCounter++;
        }

        if(globalMujocoController->controlCounter >= globalMujocoController->numControlsPerTrajectory){
            globalMujocoController->controlCounter = 0;
            globalMujocoController->resetSimFlag = true;
        }
    }

    else if(globalMujocoController->controlState == staticPos){
        double measurments[NUM_JOINTS];
        float newControls[NUM_JOINTS];

        for (int i = 0; i < NUM_JOINTS; i++) {
            measurments[i] = d->qpos[i];
        }

        for (int i = 0; i < NUM_JOINTS; i++) {
            newControls[i] = PIDController_Update(&globalMujocoController->positionPID[i], globalMujocoController->desiredJointAngles(i), measurments[i]);
        }

        for (int i = 0; i < NUM_JOINTS; i++) {
            d->ctrl[i] = newControls[i];
        }


    }
    // Linearly interpolate end effector between two 3D positions whilst countering gravity
    else if(globalMujocoController->controlState == linearInterpolation){
        m_dof testSaveControls;
        int endEffectorId = mj_name2id(model, mjOBJ_BODY, "franka_gripper");

        MatrixXd jacobian = globalMujocoController->calculateJacobian(endEffectorId);
        MatrixXd jacobian_T = jacobian.transpose();
//        cout << jacobian_T << endl;
        pose currentEndEffector = globalMujocoController->returnEndEffectorPos();

        float endEffecPosDiff[3];
        endEffecPosDiff[0] = currentEndEffector.pos.x - globalMujocoController->startEndEffectorPos.pos.x;
        endEffecPosDiff[1] = currentEndEffector.pos.y - globalMujocoController->startEndEffectorPos.pos.y;
        endEffecPosDiff[2] = currentEndEffector.pos.z - globalMujocoController->startEndEffectorPos.pos.z;

        float percentageAchieved[3];
        percentageAchieved[0] = (endEffecPosDiff[0] / globalMujocoController->diffPoseStartDesired.pos.x) * 100;
        percentageAchieved[1] = (endEffecPosDiff[1] / globalMujocoController->diffPoseStartDesired.pos.y) * 100;
        percentageAchieved[2] = (endEffecPosDiff[2] / globalMujocoController->diffPoseStartDesired.pos.z) * 100;

        Vector3d correctiveEndEffecForce;

        correctiveEndEffecForce(0) = (globalMujocoController->linearInterpolationDesiredForce(0));
        correctiveEndEffecForce(1) = (globalMujocoController->linearInterpolationDesiredForce(1));
        correctiveEndEffecForce(2) = (globalMujocoController->linearInterpolationDesiredForce(2));

//        correctiveEndEffecForce(0) = 0;
//        correctiveEndEffecForce(1) = -50;
//        correctiveEndEffecForce(2) = 0;

        // if the end effector is below correct height
        float redFactor = 100 * endEffecPosDiff[2];
        correctiveEndEffecForce(2) -= redFactor;
        //correctiveEndEffecForce(0) += redFactor;
        //cout << "corrective force z axis: " << redFactor << endl;


//        for(int i = 0; i < 2; i++){
//            float percentageDiff = (percentageAchieved[i+1] - percentageAchieved[0]);
//            if(percentageDiff > 2){
//                correctiveEndEffecForce(i+1) = globalMujocoController->linearInterpolationDesiredForce(i+1) + (5 * percentageDiff);
//            }
//            else if(percentageDiff < -2){
//                globalMujocoController->linearInterpolationDesiredForce(i+1) - (5 * percentageDiff);
//            }
//            else{
//                correctiveEndEffecForce(i+1) = globalMujocoController->linearInterpolationDesiredForce(i+1);
//            }
//
//        }

        VectorXd desiredJointTorques = jacobian_T * correctiveEndEffecForce;
        //int pointlessVariable = 0;

        for( int i = 0; i < NUM_JOINTS; i++){
            d->ctrl[i] = d->qfrc_bias[i] + desiredJointTorques(i);
            testSaveControls(i) = d->qfrc_bias[i] + desiredJointTorques(i);
        }

        saveControlCounter--;
        if(saveControlCounter <= 0){
            saveControlCounter = _mujoco_steps_per_dt;
            numControlsToSaveCounter++;
            if(numControlsToSaveCounter >= _numControls){
                saveControls(testSaveControls, true);
                if(poseAchieved(globalMujocoController->desiredEndEffectorPos, currentEndEffector)){
                    globalMujocoController->controlState = staticCalc;
                }
            }
            else{
                saveControls(testSaveControls, false);
            }
        }
    }
    else if(globalMujocoController->controlState == staticCalc){
        for( int i = 0; i < NUM_JOINTS; i++){
            d->ctrl[i] = d->qfrc_bias[i];
        }
    }
    else{

    }
}

void MujocoController::iLQRSetControlSequence(m_dof *U, int numControls){

    controlSequence = new m_dof[numControls];
    for(int i = 0; i < numControls; i++){
        controlSequence[i] = U[i];
    }
}

void initialseController(){
    mjcb_control = myController;

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

void saveControls(m_dof lastControls, bool fin){
    saveControlsFile << lastControls(0) << "," << lastControls(1) << "," << lastControls(2) << "," << lastControls(3) << "," << lastControls(4) << "," << lastControls(5) << "," << lastControls(6) << endl;
    if(fin){
        saveControlsFile.close();
    }
}

void initialiseLinearInterpolation(pose _startPose, pose _endPose, float forceMagnitude){
    globalMujocoController->startEndEffectorPos = _startPose;
    //globalMujocoController->desiredEndEffectorPos = _endPose;
    pose diff;

    diff.pos.x = globalMujocoController->desiredEndEffectorPos.pos.x - globalMujocoController->startEndEffectorPos.pos.x;
    diff.pos.y = globalMujocoController->desiredEndEffectorPos.pos.y - globalMujocoController->startEndEffectorPos.pos.y;
    diff.pos.z = globalMujocoController->desiredEndEffectorPos.pos.z - globalMujocoController->startEndEffectorPos.pos.z;
    globalMujocoController->diffPoseStartDesired = diff;
    float magnitudeDiff = sqrt(pow(diff.pos.x,2) + pow(diff.pos.y,2) + pow(diff.pos.z,2));

    diff.pos.x /= magnitudeDiff;
    diff.pos.y /= magnitudeDiff;
    diff.pos.z /= magnitudeDiff;

    globalMujocoController->linearInterpolationDesiredForce(0) = diff.pos.x * forceMagnitude;
    globalMujocoController->linearInterpolationDesiredForce(1) = diff.pos.y * forceMagnitude;
    globalMujocoController->linearInterpolationDesiredForce(2) = diff.pos.z * forceMagnitude;
    cout << "linear interpolation desired force: x:" << globalMujocoController->linearInterpolationDesiredForce(0) << " y: " << globalMujocoController->linearInterpolationDesiredForce(1) << " z: " << globalMujocoController->linearInterpolationDesiredForce(2) << endl;
    int a = 1;
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
