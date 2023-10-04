classdef RigidBodySimulation < handle
    % This class draws hawkmoth body and wings based on Morphometrics
    % objects and animates its motion based on Kinematics object 
    % wings and body can move relative to one another
    properties (Constant)
        color_body = [0.3 0.3 0.3];
        color_wing_lt = [0.5 0.9 0.1];
        color_wing_rt = [0.9 0.5 0.1];
        R_wl2m = [-1 0 0;0 -1 0;0 0 1]; % rotation matrix of left wing frame to matlab frame
        R_wr2m = [1 0 0;0 -1 0;0 0 -1]; % rotation matrix of right wing frame to matlab frame
        R_b2m = [1 0 0;0 -1 0;0 0 -1]; % rotation matrix of body frame to matlab frame
        R_e2m = [1 0 0;0 -1 0;0 0 -1]; % rotation matrix of earth frame to matlab frame
    end
    properties
        obj_fig;
        obj_body; % body object (rigid head+thorax+abdomen combined)
        obj_wing_lt; % left wing object
        obj_wing_rt; % right wing object
        obj_video; % video object to capture frames and save the simulation video 
        saveVideo;
        pos_COM_m; % position of COM in matlab frame
    end
    methods
        function obj = RigidBodySimulation(obj_morph,obj_kins,saveVideo,nameFile)
            obj.saveVideo = saveVideo;

            % initialize body and wings assuming body origin at matlab
            % origin
            temp_pos_orig_body_m = [0;0;0];

            % wing hinge positions in body frame
            temp_pos_hinge_lt_b = [obj_morph.l1*cosd(10);-obj_morph.w_b/2;-obj_morph.l1*sind(10)]; % in body attached frame (fixed)
            temp_pos_hinge_rt_b = [obj_morph.l1*cosd(10);obj_morph.w_b/2;-obj_morph.l1*sind(10)]; % in body attached frame (fixed)
            temp_pos_hinge_lt_m = temp_pos_orig_body_m + obj.R_b2m*temp_pos_hinge_lt_b;
            temp_pos_hinge_rt_m = temp_pos_orig_body_m + obj.R_b2m*temp_pos_hinge_rt_b;

            %%% initialize wing and body objects
            obj.obj_fig = figure;hold on;
            % left wing object
            obj.obj_wing_lt = BodyPart('S',"left wing",...
                [obj_morph.x_LE flip(obj_morph.x_TE)],[obj_morph.r_i flip(obj_morph.r_i)],0*[obj_morph.r_i obj_morph.r_i],...
                temp_pos_hinge_lt_m,obj.R_wl2m,[],obj_morph.l_b/20,obj.color_wing_lt,'on','off');
            % right wing object
            obj.obj_wing_rt = BodyPart('S',"right wing",...
                [obj_morph.x_LE flip(obj_morph.x_TE)],-[obj_morph.r_i flip(obj_morph.r_i)],0*[obj_morph.r_i obj_morph.r_i],...
                temp_pos_hinge_rt_m,obj.R_wr2m,[],obj_morph.l_b/20,obj.color_wing_rt,'on','off');
            % body object. Add both wing objects as its children.
            obj.obj_body = BodyPart('V',"body",...
                obj_morph.l_b/2,obj_morph.w_b/2,obj_morph.w_b/2,...
                temp_pos_orig_body_m,obj.R_b2m,[obj.obj_wing_lt obj.obj_wing_rt],obj_morph.l_b/10,obj.color_body,'on','off');
            
            % set axes 
            xlabel('x');ylabel('y');zlabel('z');
            axis equal; grid on;title(nameFile,'Interpreter','none');
            axis([-0.1 0.3 -0.1 0.1 -0.1 0.1]);

            % camera view specifying azimuthal and elevation angles
            view(-30,20);
            view(0,0);

            % set center of mass of body + wings rigid body as position vector
            % in kinematics at t = 0
            obj.pos_COM_m = obj.R_e2m*[obj_kins.x(1);obj_kins.y(1);obj_kins.z(1)];
            
            % now shift body and wings according to the combined center of mass
            obj.obj_body = obj.obj_body.partTransform([1;0;0],0,[0;0;0],-obj.obj_body.posPartCOM);
            pos_body_COM_new = (obj.pos_COM_m - obj_morph.m_w/2*obj.R_b2m*temp_pos_hinge_lt_b - obj_morph.m_w/2*obj.R_b2m*temp_pos_hinge_rt_b)/(obj_morph.m_b + obj_morph.m_w);
            obj.obj_body = obj.obj_body.partTransform([1;0;0],0,[0;0;0],pos_body_COM_new);
         
            % rotate wings relative to body with correct chi + beta 
            % (basically set up the stroke plane)
            obj.obj_wing_lt = obj.obj_wing_lt.partTransform(obj.obj_wing_lt.matPartAxes(:,2),-obj_kins.chi(1)-obj_kins.beta(1),obj.obj_wing_lt.posPartOrigin,[0;0;0]);
            obj.obj_wing_rt = obj.obj_wing_rt.partTransform(obj.obj_wing_rt.matPartAxes(:,2),-obj_kins.chi(1)-obj_kins.beta(1),obj.obj_wing_rt.posPartOrigin,[0;0;0]);
            obj.obj_wing_lt = obj.obj_wing_lt.partTransform(obj.obj_wing_lt.matPartAxes(:,1),obj_kins.beta_roll(1),obj.obj_wing_lt.posPartOrigin,[0;0;0]);
            obj.obj_wing_rt = obj.obj_wing_rt.partTransform(obj.obj_wing_rt.matPartAxes(:,1),obj_kins.beta_roll(1),obj.obj_wing_rt.posPartOrigin,[0;0;0]);
            
            % reverse translation of COM because it will be applied once
            % again in stepTransformationsForward()
            obj.obj_body = obj.obj_body.partTransform([1;0;0],0,[0;0;0],-obj.pos_COM_m);
            % apply first step transformation to body translation, body rpy and wing rpy
            obj = obj.stepTransformationsForward(obj_kins,1);

            drawnow();
            pause;
            if obj.saveVideo
                obj.obj_video = VideoWriter(strcat(nameFile,'.avi')); %// initialize the VideoWriter object
                open(obj.obj_video);
            end
        end
        function obj = runSimulation(obj,obj_kins)
            for i_N = 2:obj_kins.N
                obj = stepTransformationsReverse(obj,obj_kins,i_N-1);
                obj = stepTransformationsForward(obj,obj_kins,i_N);
                drawnow();
                if obj.saveVideo
                    curFrame = getframe(obj.obj_fig);
                    writeVideo(obj.obj_video,curFrame);
                end
            end
            if obj.saveVideo
                close(obj.obj_video);
            end
        end
        function obj = stepTransformationsForward(obj,obj_kins,i)
            % COM translation
            obj.pos_COM_m = obj.R_e2m*[obj_kins.x(i);obj_kins.y(i);obj_kins.z(i)];
            obj.obj_body = obj.obj_body.partTransform([1;0;0],0,[0;0;0],obj.pos_COM_m);
            
            % apply roll, pitch and yaw of the body 
            obj.obj_body = obj.obj_body.partTransform(obj.obj_body.matPartAxes(:,3),obj_kins.phi_b(i),obj.pos_COM_m,[0;0;0]);
            obj.obj_body = obj.obj_body.partTransform(obj.obj_body.matPartAxes(:,2),obj_kins.chi(i),obj.pos_COM_m,[0;0;0]);
            obj.obj_body = obj.obj_body.partTransform(obj.obj_body.matPartAxes(:,1),obj_kins.psi_b(i),obj.pos_COM_m,[0;0;0]);
            
            % apply roll, pitch and yaw of the wings 
            obj.obj_wing_lt = obj.obj_wing_lt.partTransform(obj.obj_wing_lt.matPartAxes(:,3),obj_kins.phi_L(i),obj.obj_wing_lt.posPartOrigin,[0;0;0]);
            obj.obj_wing_rt = obj.obj_wing_rt.partTransform(obj.obj_wing_rt.matPartAxes(:,3),obj_kins.phi_R(i),obj.obj_wing_rt.posPartOrigin,[0;0;0]);
            obj.obj_wing_lt = obj.obj_wing_lt.partTransform(obj.obj_wing_lt.matPartAxes(:,1),obj_kins.theta_L(i),obj.obj_wing_lt.posPartOrigin,[0;0;0]);
            obj.obj_wing_rt = obj.obj_wing_rt.partTransform(obj.obj_wing_rt.matPartAxes(:,1),obj_kins.theta_R(i),obj.obj_wing_rt.posPartOrigin,[0;0;0]);
            obj.obj_wing_lt = obj.obj_wing_lt.partTransform(obj.obj_wing_lt.matPartAxes(:,2),obj_kins.alpha_L(i),obj.obj_wing_lt.posPartOrigin,[0;0;0]);
            obj.obj_wing_rt = obj.obj_wing_rt.partTransform(obj.obj_wing_rt.matPartAxes(:,2),obj_kins.alpha_R(i),obj.obj_wing_rt.posPartOrigin,[0;0;0]);
        end
        function obj = stepTransformationsReverse(obj,obj_kins,i)
            obj.pos_COM_m = obj.R_e2m*[obj_kins.x(i);obj_kins.y(i);obj_kins.z(i)];
            obj.obj_wing_lt = obj.obj_wing_lt.partTransform(obj.obj_wing_lt.matPartAxes(:,2),-obj_kins.alpha_L(i),obj.obj_wing_lt.posPartOrigin,[0;0;0]);
            obj.obj_wing_rt = obj.obj_wing_rt.partTransform(obj.obj_wing_rt.matPartAxes(:,2),-obj_kins.alpha_R(i),obj.obj_wing_rt.posPartOrigin,[0;0;0]);
            obj.obj_wing_lt = obj.obj_wing_lt.partTransform(obj.obj_wing_lt.matPartAxes(:,1),-obj_kins.theta_L(i),obj.obj_wing_lt.posPartOrigin,[0;0;0]);
            obj.obj_wing_rt = obj.obj_wing_rt.partTransform(obj.obj_wing_rt.matPartAxes(:,1),-obj_kins.theta_R(i),obj.obj_wing_rt.posPartOrigin,[0;0;0]);
            obj.obj_wing_lt = obj.obj_wing_lt.partTransform(obj.obj_wing_lt.matPartAxes(:,3),-obj_kins.phi_L(i),obj.obj_wing_lt.posPartOrigin,[0;0;0]);
            obj.obj_wing_rt = obj.obj_wing_rt.partTransform(obj.obj_wing_rt.matPartAxes(:,3),-obj_kins.phi_R(i),obj.obj_wing_rt.posPartOrigin,[0;0;0]);
            obj.obj_body = obj.obj_body.partTransform(obj.obj_body.matPartAxes(:,1),-obj_kins.psi_b(i),obj.pos_COM_m,[0;0;0]);
            obj.obj_body = obj.obj_body.partTransform(obj.obj_body.matPartAxes(:,2),-obj_kins.chi(i),obj.pos_COM_m,[0;0;0]);
            obj.obj_body = obj.obj_body.partTransform(obj.obj_body.matPartAxes(:,3),-obj_kins.phi_b(i),obj.pos_COM_m,[0;0;0]);
            obj.obj_body = obj.obj_body.partTransform([1;0;0],0,[0;0;0],-obj.pos_COM_m);
        end
    end
end