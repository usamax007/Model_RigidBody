clear

%%%%%%%% set wing kinematics and morphometrics data folders
wingkinsDataDirectory = 'C:\Users\Usama\Dropbox (GaTech)\Research\Code\MatlabPlotApp\Data_WK_SpAvgd';
morphDataDirectory = 'C:\Users\Usama\Dropbox (GaTech)\Research\Code\MatlabPlotApp\Data_morph_SpAvgd';
%%%%%%%%

%%%%%%%%% select the species here
i_species = 8; 
%%%%%%%%% select species name from the list below
species_all = ["AF","AI","AL","AP","CA","EA","EI","HE","HL","PM","SO"];
% 1: AF, 2: AI, 3: AL, 4: AP, 5: CA, 6: EA, 7: EI, 8: HE, 9: HL, 10: PM, 11: SO
species_name = species_all(i_species);

%%%%%%%%% select model and simulation settings
useBodyMechanicsModel = true; % false if you want to run the blade element model alone without the coupled body mechanics. In this case, the trajectories plotted and simulated will be the actual trajectories measured in the experiment.  
runVideoSimulation = true; % if you want to watch the 3D simulation
saveVideo = true; % if you want to save the video of the 3D simulation
drawPlots = true; % if you want to draw the plots of wing kinematics, forces and torques and body trajectories
%%%%%%%%%

%%%%%%%%% set blade element strips and time steps here
n = 100; % number of blade element strips on a wing
N = 100; % number of time steps in a wingstroke
%%%%%%%%%

saveVideo = saveVideo*runVideoSimulation;
wingKinsFileName = strcat(species_name,'_00_00_0_FourierCoeffs_01.mat');
morphFileName = strcat(species_name,'_00_0_CW_BEoutput.csv');
bladeElement_bodyMechanics([n N], useBodyMechanicsModel,morphDataDirectory,morphFileName,...
    wingkinsDataDirectory,wingKinsFileName,drawPlots,true,species_name,runVideoSimulation,saveVideo);