function fout = bladeElement_bodyMechanics(BE_info, useBodyMechanics, morphDataDirectory, morphFileName,...
    wingkinsDataDirectory, wingKinsFileName,drawPlots,useTimeVarying,speciesName,playSim,saveSim)
% This function runs coupled blade element and body mechanics models 
% set useTimeVarying = true for beta, body angles and body speeds time-varying
% set useTimeVarying = false for beta, body angles and body speeds time-averaged
% useTimeVarying does not apply to the case when useBodyMechanics = true
n = BE_info(1); % number of blade element strips
N = BE_info(2); % number of sample time steps in one wing beat

disp(speciesName);
%% load morphological params
obj_morph = Morphometrics(morphDataDirectory,morphFileName,n);

%% load kinematics data
obj_kins = Kinematics(wingkinsDataDirectory,wingKinsFileName,N,useTimeVarying);

%d_wingkins = [-0.6867	12.15717 0 -5.45146	0	-0.76991	3.21895	2.47064 2.74252 6.73326 0.85976	5.69342]; % AL
%d_wingkins = [0.16116	0.07114	0	4.78189	0	-1.66564	11.12349	-0.59443	5.256	15.22908	-0.51107	8.6868]; % HE
%obj_kins.v = 0*obj_kins.v;
%obj_kins.u = 0*obj_kins.u;
%obj_kins.w = 0*obj_kins.w;
% obj_kins.chi = 0*obj_kins.chi;
% obj_kins.beta = 0*obj_kins.beta;
% obj_kins.phi_dot_L = 0*obj_kins.beta;
% obj_kins.phi_dot_R = 0*obj_kins.beta;
% obj_kins.phi_ddot_L = 0*obj_kins.beta;
% obj_kins.phi_ddot_R = 0*obj_kins.beta;
% obj_kins.theta_dot_L = 0*obj_kins.beta;
% obj_kins.theta_dot_R = 0*obj_kins.beta;
% obj_kins.theta_ddot_L = 0*obj_kins.beta;
% obj_kins.theta_ddot_R = 0*obj_kins.beta;
% obj_kins.alpha_dot_L = 0*obj_kins.beta;
% obj_kins.alpha_dot_R = 0*obj_kins.beta;
% obj_kins.alpha_ddot_L = 0*obj_kins.beta;
% obj_kins.alpha_ddot_R = 0*obj_kins.beta;

% obj_kins.beta_dot = 0*obj_kins.beta;
% obj_kins.q = 0*obj_kins.q;
% obj_kins.beta_ddot = 0*obj_kins.beta;
% obj_kins.chi = rescaleSinusoid(obj_kins.chi,mean(obj_kins.chi)+d_wingkins(2)*pi/180,max(obj_kins.chi)-min(obj_kins.chi)+d_wingkins(3)*pi/180);
% obj_kins.beta = rescaleSinusoid(obj_kins.beta,mean(obj_kins.beta)+d_wingkins(4)*pi/180,max(obj_kins.beta)-min(obj_kins.beta)+d_wingkins(5)*pi/180);
% obj_kins.beta_roll = obj_kins.beta_roll + d_wingkins(6)*pi/180;
% obj_kins.phi_L = rescaleSinusoid(obj_kins.phi_L,mean(obj_kins.phi_L)+d_wingkins(7)*pi/180,max(obj_kins.phi_L)-min(obj_kins.phi_L)+d_wingkins(8)*pi/180);
% obj_kins.phi_R = rescaleSinusoid(obj_kins.phi_R,mean(obj_kins.phi_R)+d_wingkins(7)*pi/180,max(obj_kins.phi_R)-min(obj_kins.phi_R)+d_wingkins(8)*pi/180);
% obj_kins.alpha_L = rescaleSinusoid(obj_kins.alpha_L,mean(obj_kins.alpha_L)+d_wingkins(9)*pi/180,max(obj_kins.alpha_L)-min(obj_kins.alpha_L)+d_wingkins(10)*pi/180);
% obj_kins.alpha_R = rescaleSinusoid(obj_kins.alpha_R,mean(obj_kins.alpha_R)+d_wingkins(9)*pi/180,max(obj_kins.alpha_R)-min(obj_kins.alpha_R)+d_wingkins(10)*pi/180);
% obj_kins.theta_L = rescaleSinusoid(obj_kins.theta_L,mean(obj_kins.theta_L)+d_wingkins(11)*pi/180,max(obj_kins.theta_L)-min(obj_kins.theta_L)+d_wingkins(12)*pi/180);
% obj_kins.theta_R = rescaleSinusoid(obj_kins.theta_R,mean(obj_kins.theta_R)+d_wingkins(11)*pi/180,max(obj_kins.theta_R)-min(obj_kins.theta_R)+d_wingkins(12)*pi/180);

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

%% draw plots of wing kinematics, forces and torques and body trajectories
if drawPlots
    drawPlots_bladeElement(speciesName,obj_kins,obj_morph,obj_BE_L,obj_BE_R);
end

%% Play simulation of the rigid body model
if playSim
    % initialize the rigid body simulation
    obj_sim = RigidBodySimulation(obj_morph,obj_kins,saveSim,strcat('video_',speciesName,'_test'));
    % run the simulation
    obj_sim.runSimulation(obj_kins);
end
fout = true;
