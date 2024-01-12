classdef TrimSearch
    properties
        trimmed = -1;
        tolTrim = 1e-6; 
        params_kin_min; % lower bounds of param search space
        params_kin_max; % upper bounds of param search space
        params_kin_ini; % initial condition of param search space
        params_kin_inp; % values of params from experimental data
        params_kin_trm; % trimmed values of params
        % Order of params: [f beta_roll a_chi(1)+a_beta(1) a_phi a_theta a_alpha b_phi b_theta b_alpha]
        % Assumes fixed angle between chi and beta and beta follows chi
        % variation
        listFilenames_kinBounds; % list of filenames used to determine param bounds
        n_listFilenames_kinBounds; % number of files used to determine param bounds
        objective_cost; % objective (cost) function of the search
        objective_cost_ref; % reference objective (cost) of the search
        BE_info; % [n N] blade element time steps and number of strips 
        useBodyMechanicsModel; % whether to use body mechanics as well or just run with blade element model only
        folderMorphData; % folder to load morphological data from
        filenameMorph; % name of morphological data file name
        folderWingkinsData; % folder containing averaged wing kinematics data from experiment
        folderWingkinsDataRaw; % folder containing all the raw wing kinematics data from experiments, the data used to determine bounds
        fileNameWingkins; % file name of averaged wing kinematics data from experiment
        fileNameWingkinsDummy = "kinsKinsTemp.mat"; % file name of trimmed wing kinematics data from the search. This file is created in this class.
        nameSpecies;
    end
    methods
        % constructor
        function obj = TrimSearch(numIter,costRef,BE_info, useBodyMechanics, folderMorphData, morphFileName,...
                    folderWingkinsData,folderWingkinsDataRaw, wingKinsFileName)
            obj.objective_cost_ref = costRef;
            temp_name = char(wingKinsFileName);
            obj.nameSpecies = temp_name(1:2);
            obj.BE_info = BE_info;
            obj.useBodyMechanicsModel = useBodyMechanics;
            obj.folderMorphData = folderMorphData;
            obj.filenameMorph = morphFileName;
            obj.folderWingkinsData = folderWingkinsData;
            obj.folderWingkinsDataRaw = folderWingkinsDataRaw;
            obj.fileNameWingkins = wingKinsFileName;
            obj = obj.initializeKinParams(); % initialize search bounds, initial conditions and the experimental data 
            % set up options for fminsearchbnd() function
            options = optimset('PlotFcns','optimplotfval','MaxIter',numIter,'TolFun',1e-10);
            [obj.params_kin_trm,obj.objective_cost] = fminsearchbnd(@obj.Trim_bladeElement_bodyMechanics,obj.params_kin_ini,obj.params_kin_min,obj.params_kin_max,options);
            
            if obj.objective_cost <= obj.objective_cost_ref
                obj.trimmed = true;
            end
            % plot kinematic bounds as well as initial kinematics, trimmed
            % kinematics and those from the experiment
            obj.plotKinematicBounds();

            % use trimmed kinematics to create dummy kinematics data which
            % can be input to the model
            obj.createDummyKinsData(obj.params_kin_trm);
        end
        function cost_cur = Trim_bladeElement_bodyMechanics(obj,params)
                obj.createDummyKinsData(params); % create new file with each search iterations with new params
                % run blade element and body mechanics model
                [~,F_out] = bladeElement_rigidBodyMechanics(obj.BE_info, obj.useBodyMechanicsModel,obj.folderMorphData,obj.filenameMorph,...
                    cd,obj.folderWingkinsDataRaw,obj.fileNameWingkinsDummy,false,true,obj.nameSpecies,false,false,false);
                
    

                % calculate objective (cost) function
                F_total_normd = F_out.F_total_L_normd_h + F_out.F_total_R_normd_h;
                M_total_normd = F_out.M_total_L_normd_h + F_out.M_total_R_normd_h;
                F_mean_net_normd = mean(F_total_normd,2) + [0;0;1]; % add normalized weight
                M_mean_net_normd = mean(M_total_normd,2);
                cost_cur = F_mean_net_normd'*F_mean_net_normd + M_mean_net_normd'*M_mean_net_normd;
        end
        function obj = initializeKinParams(obj)
            % Process main kinematics file
            load(strcat(obj.folderWingkinsData,'\',obj.fileNameWingkins),"f","a_chi","a_beta","beta_roll",...
                    "a_phi","b_phi","a_theta","b_theta","a_alpha","b_alpha");
            obj.params_kin_inp = [f beta_roll a_chi(1)+a_beta(1) a_phi a_theta a_alpha b_phi b_theta b_alpha];

            % Process raw kinematic files for find min and max bounds on
            % kinematic parameters
            disp('----------------------------------------------------------');
            disp('Raw kinematics files loaded to determine trim search bounds:')
            fileObjsWingKins_all = dir(obj.folderWingkinsDataRaw); % read files in the dir
            n_filesWingKins = length(fileObjsWingKins_all) - 2; % had to do this for some weird reason
            maxLengthChar = 0;
            % find maximum char length among all file names
            for i = 1:n_filesWingKins
                if(length(fileObjsWingKins_all(i+2).name)>maxLengthChar)
                    maxLengthChar = length(fileObjsWingKins_all(i+2).name);
                end
            end
            obj.listFilenames_kinBounds = char(ones(0,maxLengthChar)*'#'); % initialize char list of file names
            counterFile = 0;
            for i = 1:n_filesWingKins
                if (strcmp(obj.nameSpecies,fileObjsWingKins_all(i+2).name(1:2))) % if first two letters of the file match with the species name
                    counterFile = counterFile + 1;
                    for j = 1:length(fileObjsWingKins_all(i+2).name) % for all characters in the name
                        obj.listFilenames_kinBounds(counterFile,j) = fileObjsWingKins_all(i+2).name(j); % add char of the filename to the list
                    end
                    disp(obj.listFilenames_kinBounds(counterFile,:)); % display the filename added
                end
            end
            disp('----------------------------------------------------------');
            obj.n_listFilenames_kinBounds = counterFile;
            % order of a_params and b_params: [a_phi;a_theta;a_alpha]
            f_max = -inf;
            a_params_max = -inf*ones(3,4);
            b_params_max = -inf*ones(3,4);
            beta_roll_max = -inf;
            gamma_max = -inf; % angle between body pitch angle and stroke plane angle
            f_min = inf;
            a_params_min = inf*ones(3,4);
            b_params_min = inf*ones(3,4);
            beta_roll_min = inf;
            gamma_min = inf; % angle between body pitch angle and stroke plane angle
            for i_file = 1:obj.n_listFilenames_kinBounds
                fileName_cur = obj.listFilenames_kinBounds(i_file,:);
                load(strcat(obj.folderWingkinsDataRaw,'\',fileName_cur),"f","a_chi","a_beta","beta_roll",...
                    "a_phi","b_phi","a_theta","b_theta","a_alpha","b_alpha");
                if f_max < f
                    f_max = f;
                end
                if f_min > f
                    f_min = f;
                end
                for i_param = 1:4
                    if a_params_max(1,i_param) < a_phi(i_param)
                        a_params_max(1,i_param) = a_phi(i_param);
                    end
                    if b_params_max(1,i_param) < b_phi(i_param)
                        b_params_max(1,i_param) = b_phi(i_param);
                    end
                    if a_params_max(2,i_param) < a_theta(i_param)
                        a_params_max(2,i_param) = a_theta(i_param);
                    end
                    if b_params_max(2,i_param) < b_theta(i_param)
                        b_params_max(2,i_param) = b_theta(i_param);
                    end
                    if a_params_max(3,i_param) < a_alpha(i_param)
                        a_params_max(3,i_param) = a_alpha(i_param);
                    end
                    if b_params_max(3,i_param) < b_alpha(i_param)
                        b_params_max(3,i_param) = b_alpha(i_param);
                    end
                    if a_params_min(1,i_param) > a_phi(i_param)
                        a_params_min(1,i_param) = a_phi(i_param);
                    end
                    if b_params_min(1,i_param) > b_phi(i_param)
                        b_params_min(1,i_param) = b_phi(i_param);
                    end
                    if a_params_min(2,i_param) > a_theta(i_param)
                        a_params_min(2,i_param) = a_theta(i_param);
                    end
                    if b_params_min(2,i_param) > b_theta(i_param)
                        b_params_min(2,i_param) = b_theta(i_param);
                    end
                    if a_params_min(3,i_param) > a_alpha(i_param)
                        a_params_min(3,i_param) = a_alpha(i_param);
                    end
                    if b_params_min(3,i_param) > b_alpha(i_param)
                        b_params_min(3,i_param) = b_alpha(i_param);
                    end
                end
                if beta_roll_max < beta_roll
                    beta_roll_max = beta_roll;
                end
                if beta_roll_min > beta_roll
                    beta_roll_min = beta_roll;
                end
                if gamma_max < a_chi(1)+a_beta(1)
                    gamma_max = a_chi(1)+a_beta(1);
                end
                if gamma_min > a_chi(1)+a_beta(1)
                    gamma_min = a_chi(1)+a_beta(1);
                end
            end
            % transpose so that flatten them later using '(:)' can work properly
            a_params_max = a_params_max';
            a_params_min = a_params_min';
            b_params_max = b_params_max';
            b_params_min = b_params_min';
            obj.params_kin_max = [f_max beta_roll_max gamma_max a_params_max(:)' b_params_max(:)'];
            obj.params_kin_min = [f_min beta_roll_min gamma_min a_params_min(:)' b_params_min(:)'];
            obj.params_kin_ini = obj.params_kin_inp;
            % generate initial param values as a random point in the search
            % space
            for i_param = 1:length(obj.params_kin_min)
                obj.params_kin_ini(i_param) = obj.params_kin_min(i_param) + rand*(obj.params_kin_max(i_param)-obj.params_kin_min(i_param));
            end
        end
        function obj = plotKinematicBounds(obj)
            obj.createDummyKinsData(obj.params_kin_min);
            obj_kins_min = Kinematics(cd,obj.fileNameWingkinsDummy,obj.BE_info(2),true,false); 
            obj.createDummyKinsData(obj.params_kin_max);
            obj_kins_max = Kinematics(cd,obj.fileNameWingkinsDummy,obj.BE_info(2),true,false);     
            obj.createDummyKinsData(obj.params_kin_ini);
            obj_kins_ini = Kinematics(cd,obj.fileNameWingkinsDummy,obj.BE_info(2),true,false);
            obj.createDummyKinsData(obj.params_kin_inp);
            obj_kins_inp = Kinematics(cd,obj.fileNameWingkinsDummy,obj.BE_info(2),true,false); 
            obj.createDummyKinsData(obj.params_kin_trm);
            obj_kins_trm = Kinematics(cd,obj.fileNameWingkinsDummy,obj.BE_info(2),true,false); 
            
            % plot wing kinematics
            figure;
            % left wing kinematics
            sgtitle('Solid line: trimmed, Broken line: initial, Dotted line: averaged experimental data');
            subplot(221);hold on;
            value_trans = 0.2;
            plot(obj_kins_trm.t_norm,obj_kins_trm.phi_L*180/pi,'b','linewidth',2);
            plot(obj_kins_trm.t_norm,obj_kins_trm.alpha_L*180/pi,'r','linewidth',2);
            plot(obj_kins_trm.t_norm,obj_kins_trm.theta_L*180/pi,'k','linewidth',2);
            plot(obj_kins_ini.t_norm,obj_kins_ini.phi_L*180/pi,'b--','linewidth',2);
            plot(obj_kins_ini.t_norm,obj_kins_ini.alpha_L*180/pi,'r--','linewidth',2);
            plot(obj_kins_ini.t_norm,obj_kins_ini.theta_L*180/pi,'k--','linewidth',2);
            plot(obj_kins_inp.t_norm,obj_kins_inp.phi_L*180/pi,'b:','linewidth',2);
            plot(obj_kins_inp.t_norm,obj_kins_inp.alpha_L*180/pi,'r:','linewidth',2);
            plot(obj_kins_inp.t_norm,obj_kins_inp.theta_L*180/pi,'k:','linewidth',2);
            patch('XData',[obj_kins_min.t_norm flip(obj_kins_max.t_norm)],'YData',[obj_kins_min.phi_L*180/pi flip(obj_kins_max.phi_L*180/pi)],'FaceColor','b','FaceAlpha',value_trans,'EdgeColor','none');
            patch('XData',[obj_kins_min.t_norm flip(obj_kins_max.t_norm)],'YData',[obj_kins_min.theta_L*180/pi flip(obj_kins_max.theta_L*180/pi)],'FaceColor','k','FaceAlpha',value_trans,'EdgeColor','none');
            patch('XData',[obj_kins_min.t_norm flip(obj_kins_max.t_norm)],'YData',[obj_kins_min.alpha_L*180/pi flip(obj_kins_max.alpha_L*180/pi)],'FaceColor','r','FaceAlpha',value_trans,'EdgeColor','none');
            %
            % right wing kinematics
            title(strcat("Left wing kinematics of species ",obj.nameSpecies));
            xlabel('normalized time per wingstroke');ylabel('angle (deg)');
            axis([0 1 -90 180]);
            grid on;
            subplot(222);hold on;
            plot(obj_kins_trm.t_norm,obj_kins_trm.phi_R*180/pi,'b','linewidth',2);
            plot(obj_kins_trm.t_norm,obj_kins_trm.alpha_R*180/pi,'r','linewidth',2);
            plot(obj_kins_trm.t_norm,obj_kins_trm.theta_R*180/pi,'k','linewidth',2);
            plot(obj_kins_ini.t_norm,obj_kins_ini.phi_R*180/pi,'b--','linewidth',2);
            plot(obj_kins_ini.t_norm,obj_kins_ini.alpha_R*180/pi,'r--','linewidth',2);
            plot(obj_kins_ini.t_norm,obj_kins_ini.theta_R*180/pi,'k--','linewidth',2);
            plot(obj_kins_inp.t_norm,obj_kins_inp.phi_R*180/pi,'b:','linewidth',2);
            plot(obj_kins_inp.t_norm,obj_kins_inp.alpha_R*180/pi,'r:','linewidth',2);
            plot(obj_kins_inp.t_norm,obj_kins_inp.theta_R*180/pi,'k:','linewidth',2);
            patch('XData',[obj_kins_min.t_norm flip(obj_kins_max.t_norm)],'YData',[obj_kins_min.phi_R*180/pi flip(obj_kins_max.phi_R*180/pi)],'FaceColor','b','FaceAlpha',value_trans,'EdgeColor','none');
            patch('XData',[obj_kins_min.t_norm flip(obj_kins_max.t_norm)],'YData',[obj_kins_min.theta_R*180/pi flip(obj_kins_max.theta_R*180/pi)],'FaceColor','k','FaceAlpha',value_trans,'EdgeColor','none');
            patch('XData',[obj_kins_min.t_norm flip(obj_kins_max.t_norm)],'YData',[obj_kins_min.alpha_R*180/pi flip(obj_kins_max.alpha_R*180/pi)],'FaceColor','r','FaceAlpha',value_trans,'EdgeColor','none');
            title(strcat("Right wing kinematics of species ",obj.nameSpecies));
            xlabel('normalized time per wingstroke');ylabel('angle (deg)');
            grid on;
            legend('\phi (stroke positional angle)','\alpha (feathering angle)',...
                '\theta (stroke deviation angle)')
            axis([0 1 -90 180]);
            %
            % plot frequency ranges
            subplot(223);hold on;grid on;
            plot(obj_kins_trm.t_norm,obj_kins_trm.f*ones(1,obj_kins_trm.N),'b','linewidth',2);
            plot(obj_kins_ini.t_norm,obj_kins_ini.f*ones(1,obj_kins_ini.N),'b--','linewidth',2);
            plot(obj_kins_inp.t_norm,obj_kins_inp.f*ones(1,obj_kins_inp.N),'b:','linewidth',2);
            patch('XData',[obj_kins_min.t_norm flip(obj_kins_max.t_norm)],'YData',[obj_kins_min.f*ones(1,obj_kins_min.N) flip(obj_kins_max.f*ones(1,obj_kins_max.N))],'FaceColor','b','FaceAlpha',value_trans,'EdgeColor','none');
            xlabel('normalized time per wingstroke');ylabel('wingbeat frequency (Hz)');

            % plot beta roll and chi-beta ranges
            subplot(224);hold on;grid on;
            plot(obj_kins_trm.t_norm,obj_kins_trm.beta_roll*180/pi,'r','linewidth',2);
            plot(obj_kins_trm.t_norm,(obj_kins_trm.chi+obj_kins_trm.beta)*180/pi,'k','linewidth',2);
            plot(obj_kins_ini.t_norm,obj_kins_ini.beta_roll*180/pi,'r--','linewidth',2);
            plot(obj_kins_ini.t_norm,(obj_kins_ini.chi+obj_kins_ini.beta)*180/pi,'k--','linewidth',2);
            plot(obj_kins_inp.t_norm,obj_kins_inp.beta_roll*180/pi,'r:','linewidth',2);
            plot(obj_kins_inp.t_norm,(obj_kins_inp.chi+obj_kins_inp.beta)*180/pi,'k:','linewidth',2);
            patch('XData',[obj_kins_min.t_norm flip(obj_kins_max.t_norm)],'YData',[obj_kins_min.beta_roll*180/pi flip(obj_kins_max.beta_roll*180/pi)],'FaceColor','r','FaceAlpha',value_trans,'EdgeColor','none');
            patch('XData',[obj_kins_min.t_norm flip(obj_kins_max.t_norm)],'YData',[(obj_kins_min.chi+obj_kins_min.beta)*180/pi flip((obj_kins_max.chi+obj_kins_max.beta)*180/pi)],'FaceColor','k','FaceAlpha',value_trans,'EdgeColor','none');
            legend('\beta_{roll}','\chi+\beta');
            xlabel('normalized time per wingstroke');ylabel('angle (deg)');
            ylim([-30 150]);
        end

        function fout = createDummyKinsData(obj,params)
            % load relevant variables from the kinematics file
            load(strcat(obj.folderWingkinsData,'\',obj.fileNameWingkins),"f","a_chi","b_chi","a_beta","b_beta","beta_roll",...
                "a_phi","b_phi","a_theta","b_theta","a_alpha","b_alpha",...
                "a_u","b_u","a_v","b_v","a_w","b_w");
            f = params(1);
            beta_roll = params(2);
            a_beta(1) = params(3) - a_chi(1);
            a_beta(2:end) = -a_chi(2:end);
            b_beta = -b_chi;
            a_phi = params(4:7);
            a_theta = params(8:11);
            a_alpha = params(12:15);
            b_phi = params(16:19);
            b_theta = params(20:23);
            b_alpha = params(24:27);
            save(obj.fileNameWingkinsDummy,"f","a_chi","b_chi","a_beta","b_beta","beta_roll",...
                "a_phi","b_phi","a_theta","b_theta","a_alpha","b_alpha",...
                "a_u","b_u","a_v","b_v","a_w","b_w");
            fout = true;
        end
    end
end