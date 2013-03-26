% plots the k space reflection coefficent for the 3 layer fresnel system -
% glass, silver, vacuum.

% variables

% 633nm 
k0=9.9291803210802580537;
% indident angle (center of resonance condition) 
theta=0.54698286995313771808;
% list of thicknesses for each layer
d = [ 0 0.048068356231166835257 0 ];
% LAH79, Ag, Air
epsilon1=3.9845198023240708807+0i;
epsilon2=-14.482392074804161908+1.0945547656134573256i;
epsilon3=1+0i;
epsilon = [ epsilon1 epsilon2 epsilon3 ];
% numberof points to sample
N=1000; 
% angular spread of the region we're interested in
spread = 25.0*pi/180;

% k space
k = linspace(k0*sqrt(epsilon1)*sin(theta-spread),k0*sqrt(epsilon1)*sin(theta+spread),N);
out = nlayerfresnel(k0,k,epsilon,d);

plot(k,abs(out).^2);
