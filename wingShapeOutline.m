function [x,y_LE,y_TE] = wingShapeOutline(BE_info,morphDataDirectory,morphFileName)
n = BE_info(1);
morphData = importdata(strcat(morphDataDirectory,'\',morphFileName), ',', 1);
BEinput = morphData.data;
%wing length
data_l_w = BEinput(1,3);
%Spanwise location of each strip / blade element
%Units = meters
data_5(1,:) = BEinput(:,5) ;

%mean chord length of the each wing strip / blade element
%units = meters
data_7(1,:) = BEinput(:,7) ;

%Distance betwween the leading edge and the wing pitching axis for each strip
%Units = meters
data_8(1,:) = BEinput(:,8);

%%

[fit_c_i,S_c_i,mu_c_i] = polyfit(data_5',data_7',6);
[fit_d_w,S_d_w,mu_d_w] = polyfit(data_5',data_8',6);

%% blade element parameters
r_w_i = linspace(data_l_w/(2*n),data_l_w-data_l_w/(2*n),n)'; % distance of the midpoint of each element from the wing hinge

c_i = polyval(fit_c_i,r_w_i,S_c_i,mu_c_i)';
d_w = polyval(fit_d_w,r_w_i,S_d_w,mu_d_w)';

x = r_w_i';
y_LE = d_w;
y_TE = d_w - c_i;

end