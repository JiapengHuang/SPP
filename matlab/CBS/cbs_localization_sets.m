clc;
clear;
close all;
k = 9.9417e+06;
% number of sets of localization
NSET = 100;
%number of iterration time
NI = 100;
% distance between the input and output scatterer. The smaller distance the rougher surface is. 
dis = 0.8e-06;
field = 0;
field_abs = 0;
field_sum = 0;
rand_x = randn(1,NSET);
rand_y = randn(1,NSET);
rand_z = randn(1,NSET);
% For path 1, light is coming in the orignal point and coming out at (1,1,0.01) point. 
for j = 1:NI
    
for i=1: NSET
point_out = [dis*rand_x(i) dis*rand_y(i) 0.02e-06*rand_z(i)];
phi = linspace(-pi,pi,3600);
phase_1_in = 0;
cos_theta = 0.848;
sin_theta = 0.53;
phase_1_out = (cos_theta*point_out(3)) + (sin_theta*(cos(phi)*point_out(1) + sin(phi)*point_out(2)));
% for path 2 light is coming in the point and out at orignal.
phase_2_in = - (sin_theta*point_out(1));
phase_2_out = 0;
% show the plot
field_1 = exp(k*1i*(phase_1_in + phase_1_out));
field_2 = exp(k*1i*(phase_2_in + phase_2_out));
% here we just assume that for different sets of localization that the
% phase is same. This is also "like" the situation in the experiment since
% the phase shift is mostly determined by the mean free path. 
field = (field_1 + field_2);
field_sum = field + field_sum;
end

field_abs = field_abs + abs(field_sum);
end
max_field = max(field_abs);
field_abs = field_abs/max_field;
polar(phi,field_abs.^2);
