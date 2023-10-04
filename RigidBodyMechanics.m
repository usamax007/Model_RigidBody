classdef RigidBodyMechanics < handle
    % Aerodynamic forces acting directly on the combined body and wing
    % center of mass
    % Assuming combined body+wings as a rigid body i.e. wing movement does
    % not have inertial effects. Only aerodynamic effects.
    % The body frame is principal axis of the body. It is set by
    % rotating Matlab frame by 180 degrees about the x axis but it
    % transforms with the body transformations.
    % The earth frame is the global inertial frame. It is set by
    % rotating Matlab frame by 180 degrees about the x axis and it remains
    % fixed.
    % The RigidBodySimulation object shows that wings are moving but their
    % movement only induces aerodynamic effects and no inertial effects.
    % Center of mass of the rigid body and moment of inertia is calculated 
    % assuming both wings fully extended symmetrically and lie in the
    % body's x-y plane.
    % Moment of inertia does not change with wing movement.
    % Body drag and lift force are assumed negligible compared to
    % aerodynamic forces on the wing.
    %%%%%%%%%%%%%%%%%%%%%%%%
    % rigid body diff. equation: X_dot = A_matrix(X)*X + B_matrix(X)*F_vec(F,M,g);
    %%%%%%%%%%%%%%%%%%%%%%%%
    % A_matrix is the 12x12 nonlinear time-varying state-transition matrix
    % and is a function of the state vector X
    % B_matrix is the 12x9 nonlinear time-varying input matrix and is a
    % function of the state vector X
    % F_vec is the 9x1 input vector in which aerodynamic forces and moments, and gravitational force vectors are vertically concatenated
    properties
        X; % state vector = [x_e y_e z_e u_b v_b w_b psi_e theta_e phi_e p_b q_b r_b]
        % state vector is relative to earth frame but the subscript for each element respresents
        % the frame it was measured in (b for body and e for earth)
        % body frame is the principal axis of the body (body long axis)
        % [x y z] center of mass position
        % [u v w] center of mass velocity vector
        % [psi theta phi] roll pitch and yaw angles measured in the earth frame
        % [p q r] roll, pitch and yaw velocities measured in the body frame 
        I_b; % moment of inertia matrix in the body frame (principal axes)
        R_b2e; % rotation matrix from body frame (principal axes) to earth frame
        i_cur; % current time step index
        T; % transformation matrix for angular velocites from body rpy to earth rpy
        One; % unity matrix
        Zero; % zero matrix
        g_e; % gravity vector in the earth frame
    end
    methods
        function obj = RigidBodyMechanics(obj_morph,obj_kins)
            % constructor
            obj.i_cur = 1; % set time step index equal to 1
            obj.I_b = [obj_morph.I_xx -obj_morph.I_xy -obj_morph.I_xz;...
                -obj_morph.I_xy obj_morph.I_yy -obj_morph.I_yz;...
                -obj_morph.I_xz -obj_morph.I_yz obj_morph.I_zz]; % body frame moment of inertia. It does not vary with wing flapping.
            obj.R_b2e = rpyRotMatrix(obj_kins.psi_b(obj.i_cur),obj_kins.chi(obj.i_cur),obj_kins.phi_b(obj.i_cur));
            V_e = [obj_kins.u(obj.i_cur);obj_kins.v(obj.i_cur);obj_kins.w(obj.i_cur)]; % body velocity in earth frame
            V_b = obj.R_b2e\V_e;
            % Need to ensure kinematics p, q and r are local body p,q and r not
            % global roll, pitch and yaw rates.
            obj.X = [obj_kins.x(obj.i_cur);obj_kins.y(obj.i_cur);obj_kins.z(obj.i_cur);...
                V_b(1);V_b(2);V_b(3);...
                obj_kins.psi_b(obj.i_cur);obj_kins.chi(obj.i_cur);obj_kins.phi_b(obj.i_cur);...
                obj_kins.p(obj.i_cur);obj_kins.q(obj.i_cur);obj_kins.r(obj.i_cur)];
            obj.T = [1 sin(obj_kins.psi_b(obj.i_cur))*tan(obj_kins.chi(obj.i_cur)) cos(obj_kins.psi_b(obj.i_cur))*tan(obj_kins.chi(obj.i_cur));...
                0 cos(obj_kins.psi_b(obj.i_cur)) -sin(obj_kins.psi_b(obj.i_cur));...
                0 sin(obj_kins.psi_b(obj.i_cur))*sec(obj_kins.chi(obj.i_cur)) cos(obj_kins.psi_b(obj.i_cur))*sec(obj_kins.chi(obj.i_cur))];
            obj.One = eye(3,3);
            obj.Zero = zeros(3,3);
            obj.g_e = [0;0;9.81]; % gravity vector in earth frame
        end
        %
        function obj = stepRigidBody(obj,obj_BE_L,obj_BE_R,obj_morph,obj_kins)
            % calculate for ith time step
            Omega_b = [0 -obj.X(12) obj.X(11);...
                obj.X(12) 0 -obj.X(10);...
                -obj.X(11) obj.X(10) 0]; % skew symmteric matrix of angular speeds
            obj.R_b2e = rpyRotMatrix(obj_kins.psi_b(obj.i_cur),obj_kins.chi(obj.i_cur),obj_kins.phi_b(obj.i_cur));
            obj.T = [1 sin(obj_kins.psi_b(obj.i_cur))*tan(obj_kins.chi(obj.i_cur)) cos(obj_kins.psi_b(obj.i_cur))*tan(obj_kins.chi(obj.i_cur));...
                0 cos(obj_kins.psi_b(obj.i_cur)) -sin(obj_kins.psi_b(obj.i_cur));...
                0 sin(obj_kins.psi_b(obj.i_cur))*sec(obj_kins.chi(obj.i_cur)) cos(obj_kins.psi_b(obj.i_cur))*sec(obj_kins.chi(obj.i_cur))];
            %%
            % rigid body system diff equation: X_dot = A_matrix(X) + B_matrix(X)*F_vec(F,G,g);
            %%
            % A_matrix is the 12x12 nonlinear time-varying state-transition matrix and is a function of state variables
            % B_matrix is the 12x9 nonlinear time-varying input matrix and is a function of state variables
            % F_vec is the 9x1 input vector in which aerodynamic forces and moments, and gravitational force vectors are vertically concatenated
            A_matrix = [obj.Zero obj.R_b2e obj.Zero obj.Zero; ...
                obj.Zero -Omega_b obj.Zero obj.Zero;...
                obj.Zero obj.Zero obj.Zero obj.T;...
                obj.Zero obj.Zero obj.Zero -obj.I_b\Omega_b*obj.I_b]; % state transition matrix
            B_matrix = [obj.Zero obj.Zero obj.Zero;...
                obj.One/obj_morph.m_t obj.Zero obj.R_b2e';...
                obj.Zero obj.Zero obj.Zero;...
                obj.Zero obj.I_b\obj.One obj.Zero]; % external force and torque input matrix

            % transform torques from wing-attached to body frame
            F0_L_w = obj_BE_L.F_total_w(:,obj.i_cur);
            F0_L_b = w2b(obj_kins.phi_L(obj.i_cur),obj_kins.theta_L(obj.i_cur),...
                obj_kins.beta(obj.i_cur),obj_kins.beta_roll(obj.i_cur),...
                obj_kins.psi_b(obj.i_cur),obj_kins.chi(obj.i_cur),obj_kins.phi_b(obj.i_cur),...
                F0_L_w,'forward','L');
            F0_R_w = obj_BE_R.F_total_w(:,obj.i_cur);
            F0_R_b = w2b(obj_kins.phi_R(obj.i_cur),obj_kins.theta_R(obj.i_cur),...
                obj_kins.beta(obj.i_cur),obj_kins.beta_roll(obj.i_cur),...
                obj_kins.psi_b(obj.i_cur),obj_kins.chi(obj.i_cur),obj_kins.phi_b(obj.i_cur),...
                F0_R_w,'forward','R');
            F0_b = F0_L_b + F0_R_b;
            M0_L_w = obj_BE_L.M_total_w(:,obj.i_cur);
            M0_L_b = w2b(obj_kins.phi_L(obj.i_cur),obj_kins.theta_L(obj.i_cur),...
                obj_kins.beta(obj.i_cur),obj_kins.beta_roll(obj.i_cur),...
                obj_kins.psi_b(obj.i_cur),obj_kins.chi(obj.i_cur),obj_kins.phi_b(obj.i_cur),...
                M0_L_w,'forward','L');
            M0_R_w = obj_BE_R.M_total_w(:,obj.i_cur);
            M0_R_b = w2b(obj_kins.phi_R(obj.i_cur),obj_kins.theta_R(obj.i_cur),...
                obj_kins.beta(obj.i_cur),obj_kins.beta_roll(obj.i_cur),...
                obj_kins.psi_b(obj.i_cur),obj_kins.chi(obj.i_cur),obj_kins.phi_b(obj.i_cur),...
                M0_R_w,'forward','R');
            M0_b = M0_L_b + M0_R_b;

            % make F_vec
            F_vec = [F0_b;M0_b;obj.g_e];

            dX = obj_kins.dt*(A_matrix*obj.X + B_matrix*F_vec);
            obj.X = obj.X + dX;
            obj = updateKinematics(obj,obj_kins);
            obj.i_cur = obj.i_cur + 1;
        end

        function obj = updateKinematics(obj,obj_kins)
            % update positions
            obj_kins.x(obj.i_cur+1) = obj.X(1);
            obj_kins.y(obj.i_cur+1) = obj.X(2);
            obj_kins.z(obj.i_cur+1) = obj.X(3);
            obj_kins.psi_b(obj.i_cur+1) = obj.X(7);
            obj_kins.chi(obj.i_cur+1) = obj.X(8);
            obj_kins.phi_b(obj.i_cur+1) = obj.X(9);
            % update velocities
            obj.R_b2e = rpyRotMatrix(obj_kins.psi_b(obj.i_cur+1),obj_kins.chi(obj.i_cur+1),obj_kins.phi_b(obj.i_cur+1));
            V_b = [obj.X(4);obj.X(5);obj.X(6)];
            V_e = obj.R_b2e*V_b;
            obj_kins.u(obj.i_cur+1) = V_e(1);
            obj_kins.v(obj.i_cur+1) = V_e(2);
            obj_kins.w(obj.i_cur+1) = V_e(3);
            obj_kins.p(obj.i_cur+1) = obj.X(10);
            obj_kins.q(obj.i_cur+1) = obj.X(11);
            obj_kins.r(obj.i_cur+1) = obj.X(12);
            obj_kins.beta(obj.i_cur+1) = obj_kins.beta(obj.i_cur) - (obj_kins.chi(obj.i_cur+1) - obj_kins.chi(obj.i_cur));
            obj_kins.beta_dot(obj.i_cur+1) = (obj_kins.beta(obj.i_cur+1) - obj_kins.beta(obj.i_cur))/obj_kins.dt;
            obj_kins.beta_ddot(obj.i_cur+1) = (obj_kins.beta_dot(obj.i_cur+1) - obj_kins.beta_dot(obj.i_cur))/obj_kins.dt;
        end
    end
end