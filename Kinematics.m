classdef Kinematics < handle
    % x, y and z postions start from origin at t = 0 and evolve assuming
    % constant acceleration in the time-step dt.
    % u, v and w meausred in the earth frame (need to make sure it is 
    % consistent with how they are being extracted from experimental data).
    % It's a handle object so that the kinematics can be passed by reference
    % and overwritten in the RigidBodyMechanics model
    % In the constructor, right wing kinematics are just copied into left wing
    % kinematics from the experimental data on longitudinal and symmetric 
    % steady forward flight. But they can be modified outside once the 
    % Kinematics object is constructed because it is a handle object.
    properties
        filename; 
        N; % total time steps
        f; % wingbeat frequency
        T; % time period
        dt; % time step
        t; % time array
        t_norm; % normalized time array 
        phi_L; % left wing sweep angle
        phi_dot_L; % left wing sweep velocity
        phi_ddot_L; % left wing sweep acceleration
        theta_L; % left wing deviation angle
        theta_dot_L; % left wing deviation velocity
        theta_ddot_L; % left wing deviation acceleration
        alpha_L; % left wing pitch or feathering angle
        alpha_dot_L; % left wing pitch velocity
        alpha_ddot_L; % left wing pitch acceleration
        phi_R; % right wing sweep angle
        phi_dot_R; % right wing sweep velocity
        phi_ddot_R; % right wing sweep acceleration
        theta_R; % right wing deviation angle
        theta_dot_R; % right wing deviation velocity
        theta_ddot_R; % right wing deviation acceleration
        alpha_R; % right wing pitch or feathering angle
        alpha_dot_R; % right wing pitch velocity
        alpha_ddot_R; % right wing pitch acceleration
        psi_b; % body roll angle relative to earth frame
        chi; % body pitch angle relative to absolute horizontal
        phi_b; % body yaw angle relative to earth frame
        chi_2; % elevation of wing hinge relative to the absolute horizontal
        beta; % stroke plane angle
        beta_dot; % stroke plane angle variation rate
        beta_ddot; % stroke plane angle acceleration
        beta_roll; % tilt of the stroke plane about the stroke plane x axis
        x; % x position of body in earth frame
        y; % y position of body in earth frame
        z; % z position of body in earth frame
        u; % forward body speed in earth frame
        v; % sideslip body speed in earth frame
        w; % vertical body speed in earth frame
        p; % body frame roll rate
        q; % body frame pitch rate
        r; % body frame yaw rate
    end
    methods
        function obj = Kinematics(directory,file,n_timeSteps,useTimeVarying,dispProgress)
            obj.filename = file;
            % load relevant variables from the kinematics file
            load(strcat(directory,'\',file),"f","a_chi","b_chi","a_beta","b_beta","beta_roll",...
                "a_phi","b_phi","a_theta","b_theta","a_alpha","b_alpha",...
                "a_u","b_u","a_v","b_v","a_w","b_w");
            if dispProgress
                disp(strcat(directory,'\',file));
            end
            if(~useTimeVarying) % if time-averaged body kinematics
                % make all fourier coefficients zero for n > 0 to kill the variation
                a_chi(2:end) = 0*a_chi(2:end);
                a_beta(2:end) = 0*a_beta(2:end);
                a_u(2:end) = 0*a_u(2:end);
                a_v(2:end) = 0*a_v(2:end);
                a_w(2:end) = 0*a_w(2:end);
                b_chi = 0*b_chi;
                b_beta = 0*b_beta;
                b_u = 0*b_u;
                b_v = 0*b_v;
                b_w = 0*b_w;
            end
            obj.f = f;
            obj.N = n_timeSteps;
            obj.T = 1/f; % time period
            dt = obj.T/(n_timeSteps-1);
            t = linspace(0,obj.T,n_timeSteps); % time array
            obj.t_norm = t/obj.T;
            % all angles are defined in radians for calculations 

            %% Fourier Coefficients
            phi = a_phi(1);
            theta = a_theta(1);
            alpha = a_alpha(1);
            chi = a_chi(1);
            beta = a_beta(1);
            u = a_u(1);
            v = a_v(1);
            w = a_w(1);
            k_max = length(a_phi)-1;
            %% Compute wing and body kinematics from fitted coefficients
            for k = 1:k_max
                phi = phi + a_phi(k+1)*cos(2*pi*k*f*t) + b_phi(k+1)*sin(2*pi*k*f*t);
                alpha = alpha + a_alpha(k+1)*cos(2*pi*k*f*t) + b_alpha(k+1)*sin(2*pi*k*f*t);
                theta = theta + a_theta(k+1)*cos(2*pi*k*f*t) + b_theta(k+1)*sin(2*pi*k*f*t);
                chi = chi + a_chi(k+1)*cos(2*pi*k*f*t) + b_chi(k+1)*sin(2*pi*k*f*t);
                beta = beta + a_beta(k+1)*cos(2*pi*k*f*t) + b_beta(k+1)*sin(2*pi*k*f*t);
                u = u + a_u(k+1)*cos(2*pi*k*f*t) + b_u(k+1)*sin(2*pi*k*f*t);
                v = v + a_v(k+1)*cos(2*pi*k*f*t) + b_v(k+1)*sin(2*pi*k*f*t);
                w = w + a_w(k+1)*cos(2*pi*k*f*t) + b_w(k+1)*sin(2*pi*k*f*t);
            end
            u_dot = gradient(u)/dt;
            v_dot = gradient(v)/dt;
            w_dot = gradient(w)/dt;
            x = zeros(1,n_timeSteps);
            y = zeros(1,n_timeSteps);
            z = zeros(1,n_timeSteps);
            for i_step = 2:n_timeSteps
                x(i_step) = x(i_step-1) + u(i_step-1)*dt + 0.5*u_dot(i_step-1)*dt^2;
                y(i_step) = y(i_step-1) + v(i_step-1)*dt + 0.5*v_dot(i_step-1)*dt^2;
                z(i_step) = z(i_step-1) + w(i_step-1)*dt + 0.5*w_dot(i_step-1)*dt^2;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.t = t;
            obj.phi_L = phi;
            obj.theta_L = theta;
            obj.alpha_L = alpha;
            obj.phi_R = phi;
            obj.theta_R = theta;
            obj.alpha_R = alpha;
            obj.x = x;
            obj.y = y;
            obj.z = z;
            obj.u = u;
            obj.v = v;
            obj.w = w;
            obj.chi = chi;
            obj.psi_b = 0*chi;
            obj.phi_b = 0*chi;
            obj.chi_2 = chi + pi/6;
            
            % either tie beta to chi instead
            % obj.beta = mean(beta) + mean(chi) - chi;
            % or use actual beta
            obj.beta = beta;
            obj.beta_roll = beta_roll*ones(1,n_timeSteps); % change description if edit this
            %
            obj.phi_dot_L = gradient(phi)/dt;
            obj.phi_ddot_L = gradient(obj.phi_dot_L)/dt;
            obj.theta_dot_L = gradient(theta)/dt;
            obj.theta_ddot_L = gradient(obj.theta_dot_L)/dt;
            obj.alpha_dot_L = gradient(alpha)/dt;
            obj.alpha_ddot_L = gradient(obj.alpha_dot_L)/dt;
            %
            obj.phi_dot_R = gradient(phi)/dt;
            obj.phi_ddot_R = gradient(obj.phi_dot_R)/dt;
            obj.theta_dot_R = gradient(theta)/dt;
            obj.theta_ddot_R = gradient(obj.theta_dot_R)/dt;
            obj.alpha_dot_R = gradient(alpha)/dt;
            obj.alpha_ddot_R = gradient(obj.alpha_dot_R)/dt;
            %
            obj.beta_dot = gradient(beta)/dt;
            obj.beta_ddot = gradient(obj.beta_dot)/dt;
            obj.p = 0*obj.beta_dot;
            obj.q = -obj.beta_dot;
            obj.r = 0*obj.beta_dot;
            obj.dt = dt;
        end
    end
end
