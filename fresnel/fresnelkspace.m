% plots the fresnel coefficents for the three layer kretschmann system

% 980nm 
k0=6.414;

% indident angle (center of resonance condition) 
theta=0.54698286995313771808;

theta_lin = linspace(0.0,1.0,N);

% list of thicknesses for each layer
d = [ 0 0.048068356231166835257 0 ];

% BK7, Ag, Air
epsilon1=1.50779^2+0i;
epsilon2=-39.474505938656506 + 2.696810256318486i;
epsilon3=1+0i;
epsilon = [ epsilon1 epsilon2 epsilon3 ];

% number of points to sample
N=1000; 

% angular spread of the region we're interested in
spread = 25.0*pi/180;

k_theta = k0*sqrt(epsilon1)*sin(theta_lin);

% k space
k = linspace(k0*sqrt(epsilon1)*sin(theta-spread),k0*sqrt(epsilon1)*sin(theta+spread),N);
% specular direction
%out = nlayerfresnel(k0,k,fliplr(epsilon),fliplr(d));

% cone
out = nlayerfresnel(k0,k_theta,fliplr(epsilon),fliplr(d));
% normalize for the cone
out = out./sum(abs(out(:).^2));

plot(theta_lin,abs(out).^2);

hold on;