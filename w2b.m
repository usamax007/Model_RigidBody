%% transform from wing-attached of right and left wings to body-attached
%% principal axis. Also includes transformation to/from stroke plane only
function v_out = w2b(phi_w,theta_w,beta,beta_roll,psi_b,theta_b,phi_b,v_in,mode,wing)
% mode is 'forward' or 'inverse' transformation
% mode is also 'forwardSP' for w2sp and 'inverseSP' for sp2w
% To form R stroke plane, absolute horizontal plane is first pitched 
% beta radians about -stp_y axis and then rolled by beta_roll radians 
% about the new x axis.  
% both stroke plane frames are basically the wing-attached frames but
% before applying any wing angular movements phi_w and theta_w
% theta_b is chi body angle
% wing pitch rotation does not cause wing-attached axes to rotate
% left wing-attached frame is obtained by rotating right wing-attached
% frame by 180 degrees about the y axis. Hence additional R2 matrix
% multiplication for left wing
% The body-horizontal frame defined as a body-attached frame parallel to 
% absolute horizontal but rotates with yaw turns. Basically captures 
% heading direction and vertical force of the moth. 
N = size(v_in,2);
v_out = zeros(size(v_in));
phi_w_temp = phi_w;
theta_w_temp = theta_w;
beta_temp = beta;
beta_roll_temp = beta_roll;
psi_b_temp = psi_b;
theta_b_temp = theta_b;
phi_b_temp = phi_b;
% make changes due to differences in wing attached coordinate frames
% for each wing
if(strcmpi(wing,'R'))
    R2 = eye(3,3);
elseif(strcmpi(wing,'L'))
    R2 = rpyRotMatrix(0,pi,0);
else
    error('Invalid wing selected in w2b()');
end
for i = 1:N
    v_cur = v_in(:,i);
    theta_w_cur = theta_w_temp(i);
    phi_w_cur = phi_w_temp(i);
    beta_cur = beta_temp(i);
    beta_roll_cur = beta_roll_temp(i);
    psi_b_cur = psi_b_temp(i);
    theta_b_cur = theta_b_temp(i);
    phi_b_cur = phi_b_temp(i);

    % rotation matrix for each rotation angle
    R_theta_w = rpyRotMatrix(theta_w_cur,0,0);
    R_phi_w = rpyRotMatrix(0,0,phi_w_cur);
    R_beta = rpyRotMatrix(0,-beta_cur,0);
    R_beta_roll = rpyRotMatrix(0,beta_roll_cur,0);
    R_psi_b = rpyRotMatrix(psi_b_cur,0,0);
    R_theta_b = rpyRotMatrix(0,theta_b_cur,0);
    R_phi_b = rpyRotMatrix(0,0,phi_b_cur);

    R_w2sp = R_phi_w*R_theta_w; % wing frame to stroke plane frame
    R_sp2b = (R_psi_b'*R_theta_b'*R_phi_b')*R2*R_beta_roll*R_beta; % stroke plane frame to body-long frame
    R_sp2bh = (R_phi_b')*R2*R_beta_roll*R_beta; % stroke plane frame to body-horizontal frame

    if(strcmpi(mode,'forward'))
        R = R_sp2b*R_w2sp;
    elseif(strcmpi(mode,'inverse'))
        R = (R_sp2b*R_w2sp)';
    elseif(strcmpi(mode,'forwardSP')) % for just wing to sp frame transformation
        R = R_w2sp;
    elseif(strcmpi(mode,'inverseSP')) % for just sp to wing frame
        R = R_w2sp';
    elseif(strcmpi(mode,'forwardBH')) % for just wing to body-horizontal frame transformation
        R = R_sp2bh*R_w2sp;
    elseif(strcmpi(mode,'inverseBH')) % for just body-horizontal to wing frame
        R = (R_sp2bh*R_w2sp)';
    else
        error('Invalid mode of transformation given in function w2b()');
    end

    v_cur = R*v_cur;
    v_out(:,i) = v_cur;
end
end

