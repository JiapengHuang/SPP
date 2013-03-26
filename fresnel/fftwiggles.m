% Calculation details for "Interference of Conically Scattered Light in
% Surface Plasmon Resonance" by Aaron Webster and Frank Vollmer (Optics
% Letters 38 2).

% All distance units are in microns.  For a given propagation distance, z,
% return the detector scale (x), E_spec (Eq. 2), E_cone (Eq. 4), and a
% propagated reference Gaussian beam, E_gauss.

% This routine uses a fast Fourier transform; proper sampling is up to the
% user.

function [x,E_spec,E_cone,E_gauss] = fftwiggles(z)
	% variables

	% 633nm 
	k0=9.9291803210802580537;
	% indident angle
	theta=0.54698286995313771808;
	% list of thicknesses for each layer
	d = [ 0 0.048068356231166835257 0 ];
	% LAH79, Ag, Air
	epsilon1=3.9845198023240708807+0i;
	epsilon2=-14.482392074804161908+1.0945547656134573256i;
	epsilon3=1+0i;
	epsilon = [ epsilon1 epsilon2 epsilon3 ];
	
	% the focussed beam
	w=4;
	% size in the x dim
	wx = w/cos(theta);
	% angular spread of the focussed beam
	spread = 25.0*pi/180;
	N = 40000;
	
	% estimated sp resonance angle
	kxi = k0*sqrt(epsilon1)*sin(theta);
	% gaussian beam in k space
	gausskx = @(kx) exp(-0.25.*(kx-kxi).^2*wx.^2).*wx./sqrt(2);

	% k space
	k = linspace(k0*sqrt(epsilon1)*sin(theta-spread),k0*sqrt(epsilon1)*sin(theta+spread),N);
	
	% x and y are the actual output
	x = [0:N-1].*2.*pi./range(k);
	E_spec = ifft(1/sqrt(2*pi).*nlayerfresnel(k0,k,epsilon,d).*gausskx(k).*exp(1.0i*sqrt(k0.*k0.*epsilon1-k.*k).*z)).*abs(k(1)-k(2)).*N;
	E_cone = ifft(1/sqrt(2*pi).*nlayerfresnel(k0,k,fliplr(epsilon),fliplr(d)).*gausskx(k).*exp(1.0i*sqrt(k0.*k0.*epsilon1-k.*k).*z)).*abs(k(1)-k(2)).*N;
	E_gauss = ifft(1/sqrt(2*pi).*gausskx(k).*exp(1.0i*sqrt(k0.*k0.*epsilon1-k.*k).*z)).*abs(k(1)-k(2)).*N;
	
	E_spec = abs(E_spec).^2;
	E_cone = abs(E_cone).^2;
	% cone is way to bright, scale down
	E_cone = max(E_spec).*E_cone./max(E_cone);
	E_gauss = abs(E_gauss).^2;

	% rescale the output window (optional)
	[a,b] = showwidth(x,E_spec);
	x = x(a:b);
	E_spec = E_spec(a:b);
	E_cone = E_cone(a:b);
	E_gauss = E_gauss(a:b);

endfunction
