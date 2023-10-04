function fout = drawPlots_bladeElement(speciesName,obj_kins,obj_morph,obj_BE_L,obj_BE_R)
fout = true;
F_mean_L_w = mean(obj_BE_L.F_total_w,2);
F_mean_R_w = mean(obj_BE_R.F_total_w,2);
M_mean_L_w = mean(obj_BE_L.M_total_w,2);
M_mean_R_w = mean(obj_BE_R.M_total_w,2);

% Plot force components
y_max = 1.2*max(max(abs(obj_BE_R.F_total_w)));
figure;
subplot(131);
plot(obj_kins.t_norm,obj_BE_R.F_trans_w(1,:),'b','linewidth',2);
hold on;
plot(obj_kins.t_norm,obj_BE_R.F_rotat_w(1,:),'r','linewidth',2);
plot(obj_kins.t_norm,obj_BE_R.F_admas_w(1,:),'m','linewidth',2);
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5);
grid on;xlabel('normalized time per wingstroke');ylabel('Fx (N)');
axis([0 1 -y_max*1.2 y_max*1.2]);
hold off;
subplot(132);
title(strcat("Force components (wing frame) of species ",speciesName," at ",num2str(mean(obj_kins.u))," m/s fore-aft airspeed"));
hold on;
plot(obj_kins.t_norm,obj_BE_R.F_trans_w(2,:),'b','linewidth',2);
plot(obj_kins.t_norm,obj_BE_R.F_rotat_w(2,:),'r','linewidth',2);
plot(obj_kins.t_norm,obj_BE_R.F_admas_w(2,:),'m','linewidth',2);
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5);
grid on;xlabel('normalized time per wingstroke');ylabel('Fy (N)');
legend('Trans','Rotat', 'Admas');
axis([0 1 -y_max*1.2 y_max*1.2]);
hold off;
subplot(133);
plot(obj_kins.t_norm,obj_BE_R.F_trans_w(3,:),'b','linewidth',2);
hold on;
plot(obj_kins.t_norm,obj_BE_R.F_rotat_w(3,:),'r','linewidth',2);
plot(obj_kins.t_norm,obj_BE_R.F_admas_w(3,:),'m','linewidth',2);
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5);
grid on;xlabel('normalized time per wingstroke');ylabel('Fz (N)');
axis([0 1 -y_max*1.2 y_max*1.2]);
hold off;

% Plot Forces xyz
y_max = 1.2*max(max(abs(obj_BE_L.F_total_w)));
figure;
subplot(231);
plot(obj_kins.t_norm,obj_BE_L.F_total_w(1,:),'b','linewidth',2);
hold on;
plot(obj_kins.t_norm,F_mean_L_w(1)*ones(1,obj_kins.N),'b--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,obj_BE_R.F_total_w(1,:),'r','linewidth',2);
plot(obj_kins.t_norm,F_mean_R_w(1)*ones(1,obj_kins.N),'r--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5,'HandleVisibility','off');
grid on;legend('Left Wing','Right Wing');
xlabel('normalized time per wingstroke');ylabel('Fx (N)');
axis([0 1 -y_max*1.2 y_max*1.2]);
hold off;
subplot(232);
title(strcat("Total aerodynamic forces (wing frame) of species ",speciesName," at ",num2str(mean(obj_kins.u))," m/s fore-aft airspeed"));
hold on;
plot(obj_kins.t_norm,obj_BE_L.F_total_w(2,:),'b','linewidth',2);
plot(obj_kins.t_norm,F_mean_L_w(2)*ones(1,obj_kins.N),'b--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,obj_BE_R.F_total_w(2,:),'r','linewidth',2);
plot(obj_kins.t_norm,F_mean_R_w(2)*ones(1,obj_kins.N),'r--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5,'HandleVisibility','off');
grid on;
xlabel('normalized time per wingstroke');ylabel('Fy (N)');
axis([0 1 -y_max*1.2 y_max*1.2]);
hold off;
subplot(233);
plot(obj_kins.t_norm,obj_BE_L.F_total_w(3,:),'b','linewidth',2);
hold on;
plot(obj_kins.t_norm,F_mean_L_w(3)*ones(1,obj_kins.N),'b--','linewidth',1.5);
plot(obj_kins.t_norm,obj_BE_R.F_total_w(3,:),'r','linewidth',2);
plot(obj_kins.t_norm,F_mean_R_w(3)*ones(1,obj_kins.N),'r--','linewidth',1.5);
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5);
grid on;
xlabel('normalized time per wingstroke');ylabel('Fz (N)');
axis([0 1 -y_max*1.2 y_max*1.2]);
% Plot Moments xyz
y_max = 1.2*max(max(abs(obj_BE_L.M_total_w)));
subplot(234);
plot(obj_kins.t_norm,obj_BE_L.M_total_w(1,:),'b','linewidth',2);
hold on;
plot(obj_kins.t_norm,M_mean_L_w(1)*ones(1,obj_kins.N),'b--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,obj_BE_R.M_total_w(1,:),'r','linewidth',2);
plot(obj_kins.t_norm,M_mean_R_w(1)*ones(1,obj_kins.N),'r--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5,'HandleVisibility','off');
grid on;legend('Left Wing','Right Wing');
xlabel('normalized time per wingstroke');ylabel('Mx (Nm)');
axis([0 1 -y_max*1.2 y_max*1.2]);
hold off;
subplot(235);
title(strcat("Total aerodynamic moment (wing frame) of species ",speciesName," at ",num2str(mean(obj_kins.u))," m/s fore-aft airspeed"));
hold on;
plot(obj_kins.t_norm,obj_BE_L.M_total_w(2,:),'b','linewidth',2);
plot(obj_kins.t_norm,M_mean_L_w(2)*ones(1,obj_kins.N),'b--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,obj_BE_R.M_total_w(2,:),'r','linewidth',2);
plot(obj_kins.t_norm,M_mean_R_w(2)*ones(1,obj_kins.N),'r--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5,'HandleVisibility','off');
grid on;
xlabel('normalized time per wingstroke');ylabel('My (Nm)');
axis([0 1 -y_max*1.2 y_max*1.2]);
hold off;
subplot(236);
plot(obj_kins.t_norm,obj_BE_L.M_total_w(3,:),'b','linewidth',2);
hold on;
plot(obj_kins.t_norm,M_mean_L_w(3)*ones(1,obj_kins.N),'b--','linewidth',1.5);
plot(obj_kins.t_norm,obj_BE_R.M_total_w(3,:),'r','linewidth',2);
plot(obj_kins.t_norm,M_mean_R_w(3)*ones(1,obj_kins.N),'r--','linewidth',1.5);
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5);
grid on;
xlabel('normalized time per wingstroke');ylabel('Mz (Nm)');
axis([0 1 -y_max*1.2 y_max*1.2]);

% Plots forces and torques in body-horizontal frame
F_total_L_h = w2b(obj_kins.phi_L,obj_kins.theta_L,obj_kins.beta,obj_kins.beta_roll,obj_kins.psi_b,obj_kins.chi,obj_kins.phi_b,obj_BE_L.F_total_w,'forwardBH','L');
F_mean_L_h = mean(F_total_L_h,2);
F_total_R_h = w2b(obj_kins.phi_R,obj_kins.theta_R,obj_kins.beta,obj_kins.beta_roll,obj_kins.psi_b,obj_kins.chi,obj_kins.phi_b,obj_BE_R.F_total_w,'forwardBH','R');
F_mean_R_h = mean(F_total_R_h,2);
M_total_L_h = w2b(obj_kins.phi_L,obj_kins.theta_L,obj_kins.beta,obj_kins.beta_roll,obj_kins.psi_b,obj_kins.chi,obj_kins.phi_b,obj_BE_L.M_total_w,'forwardBH','L');
M_mean_L_h = mean(M_total_L_h,2);
M_total_R_h = w2b(obj_kins.phi_R,obj_kins.theta_R,obj_kins.beta,obj_kins.beta_roll,obj_kins.psi_b,obj_kins.chi,obj_kins.phi_b,obj_BE_R.M_total_w,'forwardBH','R');
M_mean_R_h = mean(M_total_R_h,2);
y_norm_fact = obj_morph.weight;
y_max = 1.2*max(max(abs(F_total_L_h/y_norm_fact)));
figure;
subplot(231);hold on;
plot(obj_kins.t_norm,F_total_L_h(1,:)/y_norm_fact,'b','linewidth',2);
plot(obj_kins.t_norm,F_mean_L_h(1)/y_norm_fact*ones(1,obj_kins.N),'b--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,F_total_R_h(1,:)/y_norm_fact,'r','linewidth',2);
plot(obj_kins.t_norm,F_mean_R_h(1)/y_norm_fact*ones(1,obj_kins.N),'r--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5,'HandleVisibility','off');
grid on;legend('Left Wing','Right Wing');
xlabel('normalized time per wingstroke');ylabel('Fx/weight');
axis([0 1 -y_max*1.2 y_max*1.2]);
hold off;
subplot(232);
title(strcat("Normalized aerodynamic forces (body horizontal frame) of species ",speciesName," at ",num2str(mean(obj_kins.u))," m/s fore-aft airspeed"));
hold on;
plot(obj_kins.t_norm,F_total_L_h(2,:)/y_norm_fact,'b','linewidth',2);
plot(obj_kins.t_norm,F_mean_L_h(2)/y_norm_fact*ones(1,obj_kins.N),'b--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,F_total_R_h(2,:)/y_norm_fact,'r','linewidth',2);
plot(obj_kins.t_norm,F_mean_R_h(2)/y_norm_fact*ones(1,obj_kins.N),'r--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5,'HandleVisibility','off');
grid on;
xlabel('normalized time per wingstroke');ylabel('Fy/weight');
axis([0 1 -y_max*1.2 y_max*1.2]);
hold off;
subplot(233);hold on;
plot(obj_kins.t_norm,F_total_L_h(3,:)/y_norm_fact,'b','linewidth',2);
plot(obj_kins.t_norm,F_mean_L_h(3)/y_norm_fact*ones(1,obj_kins.N),'b--','linewidth',1.5);
plot(obj_kins.t_norm,F_total_R_h(3,:)/y_norm_fact,'r','linewidth',2);
plot(obj_kins.t_norm,F_mean_R_h(3)/y_norm_fact*ones(1,obj_kins.N),'r--','linewidth',1.5);
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5);
grid on;
xlabel('normalized time per wingstroke');ylabel('Fz/weight');
axis([0 1 -y_max*1.2 y_max*1.2]);

% Plot Moments xyz
y_norm_fact = obj_morph.weight*obj_morph.r_2;
y_max = 1.2*max(max(abs(M_total_L_h/y_norm_fact)));
subplot(234);hold on;
plot(obj_kins.t_norm,M_total_L_h(1,:)/y_norm_fact,'b','linewidth',2);
plot(obj_kins.t_norm,M_mean_L_h(1)/y_norm_fact*ones(1,obj_kins.N),'b--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,M_total_R_h(1,:)/y_norm_fact,'r','linewidth',2);
plot(obj_kins.t_norm,M_mean_R_h(1)/y_norm_fact*ones(1,obj_kins.N),'r--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5,'HandleVisibility','off');
grid on;legend('Left Wing','Right Wing');
xlabel('normalized time per wingstroke');ylabel('Mx/(weight*r_2)');
axis([0 1 -y_max*1.2 y_max*1.2]);
hold off;
subplot(235);
title(strcat("Normalized aerodynamic moment (body horizontal frame) of species ",speciesName," at ",num2str(mean(obj_kins.u))," m/s fore-aft airspeed"));
hold on;
plot(obj_kins.t_norm,M_total_L_h(2,:)/y_norm_fact,'b','linewidth',2);
plot(obj_kins.t_norm,M_mean_L_h(2)/y_norm_fact*ones(1,obj_kins.N),'b--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,M_total_R_h(2,:)/y_norm_fact,'r','linewidth',2);
plot(obj_kins.t_norm,M_mean_R_h(2)/y_norm_fact*ones(1,obj_kins.N),'r--','linewidth',1.5,'HandleVisibility','off');
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5,'HandleVisibility','off');
grid on;
xlabel('normalized time per wingstroke');ylabel('My/(weight*r_2)');
axis([0 1 -y_max*1.2 y_max*1.2]);
hold off;
subplot(236);hold on;
plot(obj_kins.t_norm,M_total_L_h(3,:)/y_norm_fact,'b','linewidth',2);
plot(obj_kins.t_norm,M_mean_L_h(3)/y_norm_fact*ones(1,obj_kins.N),'b--','linewidth',1.5);
plot(obj_kins.t_norm,M_total_R_h(3,:)/y_norm_fact,'r','linewidth',2);
plot(obj_kins.t_norm,M_mean_R_h(3)/y_norm_fact*ones(1,obj_kins.N),'r--','linewidth',1.5);
plot(obj_kins.t_norm,0*obj_kins.t_norm,'k--','linewidth',0.5);
grid on;
xlabel('normalized time per wingstroke');ylabel('Mz/(weight*r_2)');
axis([0 1 -y_max*1.2 y_max*1.2]);


% plot body positions and body velocity (translational and angular)
figure;
subplot(221);hold on;
plot(obj_kins.t_norm,obj_kins.x,'linewidth',2);
plot(obj_kins.t_norm,obj_kins.y,'linewidth',2);
plot(obj_kins.t_norm,obj_kins.z,'linewidth',2);
title('Body position in earth frame');
xlabel('normalized time per wingstroke');ylabel('position (m)');
grid on;
legend('x_e','y_e', 'z_e');
hold off;
axis([0 1 -0.2 0.2]);
%
subplot(222);hold on;
plot(obj_kins.t_norm,obj_kins.psi_b*180/pi,'linewidth',2);
plot(obj_kins.t_norm,obj_kins.chi*180/pi,'linewidth',2);
plot(obj_kins.t_norm,obj_kins.phi_b*180/pi,'linewidth',2);
plot(obj_kins.t_norm,obj_kins.beta*180/pi,'linewidth',2);
title(strcat('body angles'));
xlabel('normalized time per wingstroke');ylabel('angle (deg)');
grid on;
legend('\psi_b (roll angle)', '\chi (body pitch angle)', '\phi_b (yaw angle)','\beta (stroke plane angle)');
hold off;
axis([0 1 -180 180]);
%
subplot(223);hold on;
plot(obj_kins.t_norm,obj_kins.u,'linewidth',2);
plot(obj_kins.t_norm,obj_kins.v,'linewidth',2);
plot(obj_kins.t_norm,obj_kins.w,'linewidth',2);
title('body linear velocity');
xlabel('normalized time per wingstroke');ylabel('velocity (m/s)');
grid on;
legend('u','v','w');
hold off;
axis([0 1 -4 4]);
%
subplot(224);hold on;
plot(obj_kins.t_norm,obj_kins.p,'linewidth',2);
plot(obj_kins.t_norm,obj_kins.q,'linewidth',2);
plot(obj_kins.t_norm,obj_kins.r,'linewidth',2);
title('body angular velocity');
xlabel('normalized time per wingstroke');ylabel('velocity (rad/s)');
grid on;
legend('p','q','r');
hold off;
axis([0 1 -200 200]);

% plot relative airflow speed and angle of attack for each wing strip
length = obj_morph.n;
red = [1, 0, 0];
blue = [0, 0, 1];
color_grad = [linspace(red(1),blue(1),length)', linspace(red(2),blue(2),length)', linspace(red(3),blue(3),length)'];
figure;
subplot(221);hold on;
for i_n = 1:obj_morph.n
    if i_n == 1 || i_n == obj_morph.n
        plot(obj_kins.t_norm,obj_BE_R.V_airflow_i(i_n,:),'color',color_grad(i_n,:),'linewidth',2);
    else
        plot(obj_kins.t_norm,obj_BE_R.V_airflow_i(i_n,:),'color',color_grad(i_n,:),'linewidth',2,'HandleVisibility','off');
    end
end
xlabel('normalized time per wingstroke');ylabel('right wing relative airflow speed (m/s)');
legend('proximal','distal');
grid on;
subplot(223);hold on;
for i_n = 1:obj_morph.n
    plot(obj_kins.t_norm,obj_BE_L.V_airflow_i(i_n,:),'color',color_grad(i_n,:),'linewidth',2);
end
xlabel('normalized time per wingstroke');ylabel('left wing relative airflow speed (m/s)');
grid on;
subplot(222);hold on;
for i_n = 1:obj_morph.n
    plot(obj_kins.t_norm,obj_BE_R.alpha_eff_i(i_n,:)*180/pi,'color',color_grad(i_n,:),'linewidth',2);
end
xlabel('normalized time per wingstroke');ylabel('right wing effective AoA (deg)');
ylim([-10 100]);
grid on;
subplot(224);hold on;
for i_n = 1:obj_morph.n
    plot(obj_kins.t_norm,obj_BE_L.alpha_eff_i(i_n,:)*180/pi,'color',color_grad(i_n,:),'linewidth',2);
end
xlabel('normalized time per wingstroke');ylabel('left wing effective AoA (deg)');
grid on;ylim([-10 100]);

% plot wing kinematics
figure;
subplot(121);hold on;title('left wing');
plot(obj_kins.t_norm,obj_kins.phi_L*180/pi,'linewidth',2);
plot(obj_kins.t_norm,obj_kins.alpha_L*180/pi,'r','linewidth',2);
plot(obj_kins.t_norm,(obj_kins.alpha_L-obj_kins.beta)*180/pi,'r--','linewidth',2);
plot(obj_kins.t_norm,obj_kins.theta_L*180/pi,'k','linewidth',2);
plot(obj_kins.t_norm,obj_BE_L.alpha_eff_i(floor(obj_morph.r_2/obj_morph.l_w*obj_morph.n),:)*180/pi,'g','linewidth',2);
title(strcat('Wing Kinematics of Species-',speciesName,' at: ',num2str(sqrt(mean(obj_kins.u)^2+mean(obj_kins.w)^2)),' m/s Airspeed'));
xlabel('normalized time per wingstroke');ylabel('angle (deg)');
axis([0 1 -90 180]);
grid on;
subplot(122);hold on;title('right wing');
plot(obj_kins.t_norm,obj_kins.phi_R*180/pi,'linewidth',2);
plot(obj_kins.t_norm,obj_kins.alpha_R*180/pi,'r','linewidth',2);
plot(obj_kins.t_norm,(obj_kins.alpha_R-obj_kins.beta)*180/pi,'r--','linewidth',2);
plot(obj_kins.t_norm,obj_kins.theta_R*180/pi,'k','linewidth',2);
plot(obj_kins.t_norm,obj_BE_R.alpha_eff_i(floor(obj_morph.r_2/obj_morph.l_w*obj_morph.n),:)*180/pi,'g','linewidth',2);
title(strcat('Wing Kinematics of Species-',speciesName,' at: ',num2str(sqrt(mean(obj_kins.u)^2+mean(obj_kins.w)^2)),' m/s Airspeed'));
xlabel('normalized time per wingstroke');ylabel('angle (deg)');
grid on;
legend('\phi (stroke positional angle)','\alpha (feathering angle)',...
    '\alpha - \beta (chord angle in global frame)','\theta (stroke deviation angle)',strcat('\alpha effective (airspeed =',num2str(sqrt(mean(obj_kins.u)^2+mean(obj_kins.w)^2)),' m/s)'));
axis([0 1 -90 180]);
end
