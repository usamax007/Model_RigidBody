function [success,forces_h] = bladeElement_rigidBodyMechanics(BE_info, useBodyMechanics, morphDataDirectory, morphFileName,...
    wingkinsDataDirectory,wingkinsRawDataDirectory,wingKinsFileName,drawPlots,useTimeVarying,speciesName,runTrimSearch,playSim,saveSim)
% This function runs coupled blade element and body mechanics models 
% set useTimeVarying = true for beta, body angles and body speeds time-varying
% set useTimeVarying = false for beta, body angles and body speeds time-averaged
% useTimeVarying does not apply to the case when useBodyMechanics = true

n = BE_info(1); % number of blade element strips
N = BE_info(2); % number of sample time steps in one wing beat

if runTrimSearch
    obj_trimSearch = TrimSearch(1000,1e-5,BE_info, useBodyMechanics, morphDataDirectory, morphFileName,...
                        wingkinsDataDirectory,wingkinsRawDataDirectory, wingKinsFileName);
    wingkinsDataDirectory = cd;
    wingKinsFileName = obj_trimSearch.fileNameWingkinsDummy;
end

if drawPlots
    disp(speciesName);
end
%% load morphological params
obj_morph = Morphometrics(morphDataDirectory,morphFileName,n,drawPlots);

%% load kinematics data
obj_kins = Kinematics(wingkinsDataDirectory,wingKinsFileName,N,useTimeVarying,drawPlots);
obj_kins.v = 0*obj_kins.v;
%obj_kins = suppressWingMotion(obj_kins);


%% initialize BladeElement models for left and right wings
obj_BE_L = BladeElement(N,n,'L');
obj_BE_R = BladeElement(N,n,'R');

%% initialize the rigid body model (just the model, not the simulation)
if useBodyMechanics
    obj_RB = RigidBodyMechanics(obj_morph,obj_kins);
end

%% Run the first step starting at t = 0
obj_BE_L = obj_BE_L.stepBladeElement(1,obj_morph,obj_kins);
obj_BE_R = obj_BE_R.stepBladeElement(1,obj_morph,obj_kins);
if useBodyMechanics
    obj_RB = obj_RB.stepRigidBody(obj_BE_L,obj_BE_R,obj_morph,obj_kins);
end

%% Run the remaining time steps up to N-1 steps because the rigid body mechanics model calculates quantities for the Nth step
for i_N = 2:N-1
    obj_BE_L = obj_BE_L.stepBladeElement(i_N,obj_morph,obj_kins);
    obj_BE_R = obj_BE_R.stepBladeElement(i_N,obj_morph,obj_kins);
    if useBodyMechanics
        obj_RB = obj_RB.stepRigidBody(obj_BE_L,obj_BE_R,obj_morph,obj_kins);
    end
end

forces_h.F_total_L_normd_h = w2b(obj_kins.phi_L,obj_kins.theta_L,obj_kins.beta,obj_kins.beta_roll,obj_kins.psi_b,obj_kins.chi,obj_kins.phi_b,obj_BE_L.F_total_w/obj_morph.weight,'forwardBH','L');
forces_h.M_total_L_normd_h = w2b(obj_kins.phi_L,obj_kins.theta_L,obj_kins.beta,obj_kins.beta_roll,obj_kins.psi_b,obj_kins.chi,obj_kins.phi_b,obj_BE_L.M_total_w/(obj_morph.weight*obj_morph.r_2),'forwardBH','L');
forces_h.F_total_R_normd_h = w2b(obj_kins.phi_R,obj_kins.theta_R,obj_kins.beta,obj_kins.beta_roll,obj_kins.psi_b,obj_kins.chi,obj_kins.phi_b,obj_BE_R.F_total_w/obj_morph.weight,'forwardBH','R');
forces_h.M_total_R_normd_h = w2b(obj_kins.phi_R,obj_kins.theta_R,obj_kins.beta,obj_kins.beta_roll,obj_kins.psi_b,obj_kins.chi,obj_kins.phi_b,obj_BE_R.M_total_w/(obj_morph.weight*obj_morph.r_2),'forwardBH','R');
            

%% draw plots of wing kinematics, forces and torques and body trajectories
if drawPlots
    drawPlots_bladeElement(speciesName,obj_kins,obj_morph,obj_BE_L,obj_BE_R,forces_h);
end

%% Play simulation of the rigid body model
if playSim
    % initialize the rigid body simulation
    obj_sim = RigidBodySimulation(obj_morph,obj_kins,saveSim,strcat('video_',speciesName,'_test'));
    % run the simulation
    obj_sim.runSimulation(obj_kins);
end
success = true;
