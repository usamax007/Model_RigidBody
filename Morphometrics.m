classdef Morphometrics
    % all morphometrics are output in SI units
    properties 
        filename;
        n; % number of chordwise strips the wing is divided into (blade elements)
        m_t; % total mass
        m_w; % mass of both wings
        m_b; % mass of body (head + thorax + abdomen + legs)
        I_xx; % elements of moment of total inertia matrix in the body-long frame
        I_yy;
        I_zz;
        I_xy;
        I_yz;
        I_xz;
        l_b; % body length
        w_b; % body width
        l_w; % wing length
        l1; % distance of wing hinge from center of mass
        r_2; % radius of the second moment of area of the wing (meters)
        r_i; % distance of the midpoint of each element from the wing hinge
        l_dr; % array of width of each strip
        c_i; % array of chord length of each strip
        d_w; % array of positions of leading edge on each strip
        x_LE; % array of positions of leading edge on each strip
        x_TE; % array of positions of trailing edge on each strip
        eps_i; % array of positions of half chord on each strip
        eps_i_2; % array of positions of quarter chord on each strip
        c_mean; % mean chord length
        S_w; % wing area
        weight; % total body and wing weight
        wingCoM_w; % coordinates of wing center of mass in wing frame
    end
    methods
        function obj = Morphometrics(directory,file,nStrips,dispProgress)
            % import morphological data
            obj.n = nStrips;
            obj.filename = file;
            morphData = importdata(strcat(directory,'\',file), ',', 1);
            if dispProgress
                disp(strcat(directory,'\',file))
            end
            BEinput = morphData.data;
            %%
            % Body and wing morphology
            data_m_t = BEinput(1,1); % total mass of body plus wings (grams)
            data_m_w = BEinput(1,16); % right wing mass estimate (grams)
            data_Ixx = BEinput(1,17); % Ixx estimate (grams meter^2)
            data_Iyy = BEinput(1,18); % Iyy estimate (grams meter^2)
            data_Izz = BEinput(1,19); % Izz estimate (grams meter^2)
            data_Ixz = BEinput(1,20); % Ixz estimate (grams meter^2)
            data_Ixy = 0; % Ixy estimate (grams meter^2)
            data_Iyz = 0; % Iyz estimate (grams meter^2)
            data_l_b = BEinput(1,2); % body length (meters)
            data_l_w = BEinput(1,3); % wing length (meters)
            data_S_w = BEinput(1,12); % wing area (meters^2)
            data_r2 = BEinput(1,13); % nondim radius of the second moment of area of the wing 
            data_l1 = BEinput(1,9); % distance between the wing base pivot/pitching axis and the center of mass (COM)
            
            % Wing strip morphology
            data_l_lp = BEinput(:,8); % distance between the leading edge and wing-pitching axis (meters)
            data_r_all(1,:) = BEinput(:,5); % spanwise location of each strip / blade element (meters)
            data_c_mean_all(1,:) = BEinput(:,7); % mean chord length of the each wing strip / blade element (meters)
            data_l_lp_all(1,:) = BEinput(:,8); % distance betwween the leading edge and the wing pitching axis for each strip (meters)
            data_l_lp_all_nondim(1,:) = BEinput(:,10); % nondimensional value representing the distance between the leading edge and the wing's pitching axis.
            data_l_hp_all(1,:) = BEinput(:,11); % distance betweeen half-chord and pitching axis (meters) 
            data_l_qp_all(1,:) = BEinput(:,22); % Distance betweeen quarter-chord and pitching axis (meters)
            
            %% fitting wing shape from morph parameters
            [fit_c_i,S_c_i,mu_c_i] = polyfit(data_r_all',data_c_mean_all',6);
            [fit_d_w,S_d_w,mu_d_w] = polyfit(data_r_all',data_l_lp_all',6);
            [fit_eps_i,S_eps_i,mu_eps_i] = polyfit(data_r_all',data_l_hp_all',6);
            [fit_eps_i_2,S_eps_i_2,mu_eps_i_2] = polyfit(data_r_all',data_l_qp_all',6);
            
            %% morphological params
            obj.m_t = data_m_t*1e-3; % kg
            obj.m_w = 2*data_m_w*1e-3; % total mass of both wings
            obj.m_b = obj.m_t - obj.m_w;
            obj.I_xx = data_Ixx*1e-3;
            obj.I_yy = data_Iyy*1e-3;
            obj.I_zz = data_Izz*1e-3;
            obj.I_xz = data_Ixz*1e-3;
            obj.I_xy = data_Ixy*1e-3;
            obj.I_yz = data_Iyz*1e-3;
            obj.l_b = data_l_b;
            obj.w_b = sqrt(6*obj.m_b/(pi*obj.l_b*1000)); % assuming density of water
            obj.l_w = data_l_w;
            obj.l1 = data_l1;
            obj.r_2 = data_r2*obj.l_w;
            obj.weight = obj.m_t*9.81;
            
            %% blade element parameters
            obj.r_i = linspace(obj.l_w/(2*obj.n),obj.l_w-obj.l_w/(2*obj.n),obj.n); 
            obj.l_dr = obj.l_w/obj.n; 
            obj.c_i = polyval(fit_c_i,obj.r_i,S_c_i,mu_c_i);
            obj.d_w = polyval(fit_d_w,obj.r_i,S_d_w,mu_d_w);
            obj.eps_i = polyval(fit_eps_i,obj.r_i,S_eps_i,mu_eps_i);
            obj.eps_i_2 = polyval(fit_eps_i_2,obj.r_i,S_eps_i_2,mu_eps_i_2);
            
            obj.c_mean = mean(obj.c_i);
            obj.S_w = data_S_w;
            
            % calculate centroid of a wing in wing coordinate frame
            obj.wingCoM_w = sum([-obj.c_i.*obj.eps_i;obj.c_i.*obj.r_i;zeros(size(obj.c_i))],2)/sum(obj.c_i);
            
            obj.x_LE = obj.d_w;
            obj.x_TE = obj.d_w - obj.c_i;
            
            % shift wing hinge to the proximal tip of the leading edge
            %obj.x_TE = obj.x_TE - obj.x_LE(1);
            %obj.x_LE = obj.x_LE - obj.x_LE(1);
        end
    end
end