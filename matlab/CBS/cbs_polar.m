clc;
clear;
close all;
k = 9.9417e+06;
% For path 1, light is coming in the orignal point and coming out at (1,1,0.01) point. 
point_out = [-1.0e-06 -4.0e-06 0.01e-06];
phi = linspace(-pi,pi,3600);
phase_1_in = 0;
cos_theta = 0.848;
sin_theta = 0.53;
phase_1_out = (cos_theta*point_out(3)) + (sin_theta*(cos(phi)*point_out(1) + sin(phi)*point_out(2)));

% for path 2, light is coming in the point and out at orignal.
phase_2_in = (sin_theta*point_out(1));
phase_2_out = 0;

% show the plot
field_1 = exp(k*1i*(phase_1_in + phase_1_out));
field_2 = exp(k*1i*(phase_2_in + phase_2_out));
field = field_1 + field_2;
field_abs = abs(field);
polar(phi,field_abs);