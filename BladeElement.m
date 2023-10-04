classdef BladeElement < handle
    % class for calculating aerodynamic forces on each wing 
    % forces are calculated in the wing-attached frames
    % The right wing frame is obtained by rotating matlab coordinate frame by
    % 180 degrees about its x axis
    % The left wing frame is obtained by rotating right wing frame by 180
    % degrees about its y axis
    % Wing frames are chosen to ensure that the same wing kinematic signs apply
    % to both left and right wings for symmetric movements
    properties (Constant)
        rho_air = 1.184; % density of air at 25 C (units: kg/m^3)
        mu_air = 1.849e-5; % dynamic viscosity of air at 25 C (units: kg/m/s)
    end
    properties
        nElements; % number (n) of blade element strips
        nTimeSteps; % number (N) of time steps in the kinematic waveforms
        wing_RL; % char 'L' or 'R' 
        V_airflow_i % nxN matrix of airflow speeds of each strip (rows) at each time step (columns) 
        alpha_eff_i; % nxN matrix of effective AoA of each strip (rows) at each time step (columns) 
        
        % subscript w represents measured in wing frame
        % subscript b represents measured in body frame
        F_trans_w; % 3xN matrix of total translational aerodynamic force vector 
        F_transDrag_w; % 3xN matrix of translational drag vector 
        F_transLift_w; % 3xN matrix of translational lift vector 
        F_rotat_w; % 3xN matrix of total rotational aerodynamic force vector 
        F_admas_w; % 3xN matrix of total added-mass aerodynamic force vector 
        F_total_w; % 3xN matrix of total aerodynamic force vector 
        M_trans_w; % 3xN matrix of total translational aerodynamic moment vector 
        M_transDrag_w; % 3xN matrix of translational drag moment vector 
        M_transLift_w; % 3xN matrix of translational lift moment vector 
        M_rotat_w; % 3xN matrix of rotational aerodynamic moment vector 
        M_admas_w; % 3xN matrix of added-mass aerodynamic moment vector 
        M_total_w; % 3xN matrix of total aerodynamic moment vector 
    end
    methods
        function obj = BladeElement(N,n,wing) 
            % constructor
            % N = number of time steps in the kinematic waveforms
            % n = number of blade element strips on a wing
            % wing is entered as a char 'L' for left wing or 'R' for right wing 
            obj.nTimeSteps = N;
            obj.nElements = n;
            obj.wing_RL = wing;
            obj.V_airflow_i = zeros(n,N);
            obj.alpha_eff_i = zeros(n,N);
            obj.F_trans_w = zeros(3,N);
            obj.F_transDrag_w = zeros(3,N);
            obj.F_transLift_w = zeros(3,N);
            obj.F_rotat_w = zeros(3,N);
            obj.F_admas_w = zeros(3,N);
            obj.F_total_w = zeros(3,N);
            obj.M_trans_w = zeros(3,N);
            obj.M_transDrag_w = zeros(3,N);
            obj.M_transLift_w = zeros(3,N);
            obj.M_rotat_w = zeros(3,N);
            obj.M_admas_w = zeros(3,N);
            obj.M_total_w = zeros(3,N);
        end
        function obj = stepBladeElement(obj,i,obj_morph,obj_kins)
            %% Calculate body velocity in the wing frame
            if strcmpi(obj.wing_RL,'R') % if right wing
                phi = obj_kins.phi_R(i);
                phi_dot = obj_kins.phi_dot_R(i);
                phi_ddot = obj_kins.phi_ddot_R(i);
                theta = obj_kins.theta_R(i);
                theta_dot = obj_kins.theta_dot_R(i);
                theta_ddot = obj_kins.theta_ddot_R(i);
                alpha = obj_kins.alpha_R(i);
                alpha_dot = obj_kins.alpha_dot_R(i);
                alpha_ddot = obj_kins.alpha_ddot_R(i);
                b_hat_wing = [cos(alpha);0;-sin(alpha)]; % wing trailing to leading edge vector
                n_hat_wing = [-sin(alpha);0;-cos(alpha)]; % wing surface normal vector
                r_i_vect_w = [0*obj_morph.r_i;obj_morph.r_i;0*obj_morph.r_i]; % 3xn matrix of coordinates of wing strip centers 
            elseif strcmpi(obj.wing_RL,'L') % if left wing
                phi = obj_kins.phi_L(i);
                phi_dot = obj_kins.phi_dot_L(i);
                phi_ddot = obj_kins.phi_ddot_L(i);
                theta = obj_kins.theta_L(i);
                theta_dot = obj_kins.theta_dot_L(i);
                theta_ddot = obj_kins.theta_ddot_L(i);
                alpha = obj_kins.alpha_L(i);
                alpha_dot = obj_kins.alpha_dot_L(i);
                alpha_ddot = obj_kins.alpha_ddot_L(i);
                b_hat_wing = [-cos(alpha);0;sin(alpha)];
                n_hat_wing = [sin(alpha);0;cos(alpha)];
                r_i_vect_w = [0*obj_morph.r_i;-obj_morph.r_i;0*obj_morph.r_i]; % 3xn matrix of coordinates of wing strip centers 
            else
                error('invalid wing_RL')
            end
            R_b2e = rpyRotMatrix(obj_kins.psi_b(i),obj_kins.chi(i),obj_kins.phi_b(i));
            R_e2b = R_b2e';
            V_body_e = [obj_kins.u(i);obj_kins.v(i);obj_kins.w(i)]; % body translational velocity in earth frame
            V_body_b = R_e2b*V_body_e; % body translational velocity vector
            omega_body_b = [obj_kins.p(i);obj_kins.q(i);obj_kins.r(i)]; % body rotational velocity vector
            V_body_w = w2b(phi,theta,obj_kins.beta(i),obj_kins.beta_roll(i),obj_kins.psi_b(i),obj_kins.chi(i),obj_kins.phi_b(i),V_body_b,'inverse',obj.wing_RL); % body translational velocity in wing frame
            omega_body_w = w2b(phi,theta,obj_kins.beta(i),obj_kins.beta_roll(i),obj_kins.psi_b(i),obj_kins.chi(i),obj_kins.phi_b(i),omega_body_b,'inverse',obj.wing_RL); % body rotational velocity in wing frame
           
            %% Calculate relative airflow velocity in the wing frame
            omega_wing_sp = [0; 0; phi_dot] + ...
                    w2b(phi,0,0,0,0,0,0,[theta_dot; 0; 0],'forwardSP',obj.wing_RL); % wing rotational velocity in stroke plane frame
            omega_wing_w = w2b(phi,theta,0,0,0,0,0,omega_wing_sp,'inverseSP',obj.wing_RL); % stroke plane to wing frame
            l1_b = obj_morph.l1*[cos(obj_kins.chi_2(i)-obj_kins.chi(i));0*obj_kins.chi_2(i);-sin(obj_kins.chi_2(i)-obj_kins.chi(i))]; % l1 vector in body frame
            l1_w = w2b(phi,theta,obj_kins.beta(i),obj_kins.beta_roll(i),obj_kins.psi_b(i),obj_kins.chi(i),obj_kins.phi_b(i),l1_b,'inverse',obj.wing_RL); % l1 vector in wing frame
            FD_hat_i_w = zeros(3,obj.nElements); % drag direction in wing frame
            FL_hat_i_w = zeros(3,obj.nElements); % lift direction in wing frame
            alpha_eff_res_i = zeros(1,obj.nElements); % effective angle of attack bound between 0 and 90 degs
            for j = 1:obj.nElements
                r_i_vect_w_cur = r_i_vect_w(:,j);
                V_airflow_wing_w_cur = -cross(omega_wing_w,r_i_vect_w_cur); % relative airflow just due to wing motion (v = omega x r)
                V_total_w_cur = V_airflow_wing_w_cur - (V_body_w + cross(omega_body_w,(l1_w + r_i_vect_w_cur))); % add airflow due to body's translational and rotational velocities
                obj.V_airflow_i(j,i) = norm(V_total_w_cur); % total airflow speed
                FD_hat = normalizeVect(V_total_w_cur); % Drag unit vector is in direction of the airflow
                FL_hat = (FD_hat'*n_hat_wing)*cross(cross(FD_hat,n_hat_wing),FD_hat); % lift vector is perpendicular to drag but it's direction is a complicated function of incoming airflow vector and wing surface normal vector
                FL_hat = normalizeVect(FL_hat); % normalize just in case
                FD_hat_i_w(:,j) = FD_hat;
                FL_hat_i_w(:,j) = FL_hat;
                obj.alpha_eff_i(j,i) = acos(-FD_hat'*b_hat_wing); % effective angle of attack (AoA)
                if(obj.alpha_eff_i(j,i) < pi/2) % bound the AoA between 0 and 90
                    alpha_eff_res_i(j) = obj.alpha_eff_i(j,i);
                else
                    alpha_eff_res_i(j) = (pi - obj.alpha_eff_i(j,i));
                end
            end

            %% Calculate aerodynamic forces

            % Aerodynamic coefficients from Han 2016. The ones from Kim 2015 were
            % behaving in a weird way. I seriously think there is some problem with
            % those equations.
            CL = 1.5520*sind(alpha_eff_res_i*180/pi).*cosd(alpha_eff_res_i*180/pi) + 1.725*((sind(alpha_eff_res_i*180/pi)).^2).*cosd(alpha_eff_res_i*180/pi);
            CD = 0.0596*sind(alpha_eff_res_i*180/pi).*cosd(alpha_eff_res_i*180/pi) + 3.598*(sind(alpha_eff_res_i*180/pi)).^3;
            CR = pi*(0.75 - obj_morph.d_w./obj_morph.c_i); % d_w is length from the leading edge to the wing pitching axis

            % CL = 1.07597*CL;
            % CD = 0.66061*CD;

            % CL = 1.74511*CL;
            % CD = 1.39978*CD;
            
            %CL = 1.3489;
            %CD = 0.73273;

            % lift and drag components of translational aerodynamic force
            % magnitudes for each blade element strip at ith time instance
            F_transLift_i = 1/2*obj.rho_air*(obj_morph.l_dr)*obj_morph.c_i.*CL.*(obj.V_airflow_i(:,i).^2)';
            F_transDrag_i = 1/2*obj.rho_air*(obj_morph.l_dr)*obj_morph.c_i.*CD.*(obj.V_airflow_i(:,i).^2)';
            
            % rotational and added-mass aerodynamic force magintudes for
            % all blade element strip at ith time instance
            F_rotat_i = obj.rho_air*(obj_morph.l_dr)*(obj_morph.c_i.^2).*CR.*(obj.V_airflow_i(:,i)')*(alpha_dot - obj_kins.beta_dot(i));
            F_admas_i = (pi/4*obj.rho_air*((phi_ddot*sin(alpha) + phi_dot*(alpha_dot - obj_kins.beta_dot(i))*cos(alpha))*obj_morph.r_i.*(obj_morph.c_i.^2)*obj_morph.l_dr ...
                + 1/4*obj_morph.l_dr*(alpha_ddot-obj_kins.beta_ddot(i))*(obj_morph.c_i.^3)));
            
            F_trans_i_w = zeros(3,obj.nElements);
            F_transDrag_i_w = zeros(3,obj.nElements);
            F_transLift_i_w = zeros(3,obj.nElements);
            M_trans_i_w = zeros(3,obj.nElements);
            M_transDrag_i_w = zeros(3,obj.nElements);
            M_transLift_i_w = zeros(3,obj.nElements);
            F_rotat_i_w = zeros(3,obj.nElements);
            M_rotat_i_w = zeros(3,obj.nElements);
            F_admas_i_w = zeros(3,obj.nElements);
            M_admas_i_w = zeros(3,obj.nElements);
            for j = 1:obj.nElements
                % convert force magnitudes to vectors
                if strcmpi(obj.wing_RL,'R')
                    F_rotat_i_w(:,j) = F_rotat_i(j)*[-sin(alpha); 0; -cos(alpha)];
                    F_admas_i_w(:,j) = F_admas_i(j)*[sin(alpha); 0; cos(alpha)];
                elseif strcmpi(obj.wing_RL,'L')
                    F_rotat_i_w(:,j) = F_rotat_i(j)*[sin(alpha); 0; cos(alpha)];
                    F_admas_i_w(:,j) = F_admas_i(j)*[-sin(alpha); 0; -cos(alpha)];
                end
                FD_hat = FD_hat_i_w(:,j);
                FL_hat = FL_hat_i_w(:,j);
                F_rotat_w_temp = F_rotat_i_w(:,j);
                F_admas_w_temp = F_admas_i_w(:,j);
                
                % convert force magnitudes to vectors
                F_transDrag_temp = F_transDrag_i(j).*FD_hat;
                F_transLift_temp = F_transLift_i(j).*FL_hat;
                
                % total translational aerodynamic force
                F_trans_temp = F_transDrag_temp + F_transLift_temp;
                F_transDrag_i_w(:,j) = F_transDrag_temp;
                F_transLift_i_w(:,j) = F_transLift_temp;
                F_trans_i_w(:,j) = F_trans_temp;

                % moments generated by each force about the body center of mass
                % translational force assumed to act at the quarter chord point from the leading edge
                % rotational and added-mass forces assumed to act at the half-chord point from the leading edge 
                moment_arm_q = l1_w + r_i_vect_w(:,j) + b_hat_wing*obj_morph.eps_i_2(j); % moment arm to quarter chord
                moment_arm_h = l1_w + r_i_vect_w(:,j) + b_hat_wing*obj_morph.eps_i(j); % moment arm to half chord
                M_trans_i_w(:,j) = cross(moment_arm_q,F_trans_temp);
                M_transDrag_i_w(:,j) = cross(moment_arm_q,F_transDrag_temp);
                M_transLift_i_w(:,j) = cross(moment_arm_q,F_transLift_temp);
                M_rotat_i_w(:,j) = cross(moment_arm_h,F_rotat_w_temp);
                M_admas_i_w(:,j) = cross(moment_arm_h,F_admas_w_temp);
            end

            % sum over all the blade element strips from 1 to n
            obj.F_transDrag_w(:,i) = sum(F_transDrag_i_w,2);
            obj.F_transLift_w(:,i) = sum(F_transLift_i_w,2);
            obj.F_trans_w(:,i) = sum(F_trans_i_w,2);
            obj.M_transDrag_w(:,i) = sum(M_transDrag_i_w,2);
            obj.M_transLift_w(:,i) = sum(M_transLift_i_w,2);
            obj.M_trans_w(:,i) = sum(M_trans_i_w,2);
            obj.F_rotat_w(:,i) = sum(F_rotat_i_w,2);
            obj.M_rotat_w(:,i) = sum(M_rotat_i_w,2);
            obj.F_admas_w(:,i) = sum(F_admas_i_w,2);
            obj.M_admas_w(:,i) = sum(M_admas_i_w,2);

            % total aerodynamic force and moments
            obj.F_total_w(:,i) = obj.F_trans_w(:,i) + obj.F_rotat_w(:,i) + obj.F_admas_w(:,i);
            obj.M_total_w(:,i) = obj.M_trans_w(:,i) + obj.M_rotat_w(:,i) + obj.M_admas_w(:,i);
        end
    end
end