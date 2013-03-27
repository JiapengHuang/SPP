% returns the complex n layer Fresnel reflectivity for TM polarization
function out = nlayerfresnel(k0,kx,epsilon,d)
	kz = @(kx,epsilon) sqrt(k0.^2.*epsilon-kx.^2);
	% two layer Fresnel
	r12 = @(kx,epsilon) ...
				(epsilon(2).*kz(kx,epsilon(1))-epsilon(1).*kz(kx,epsilon(2)))./ ...
				(epsilon(2).*kz(kx,epsilon(1))+epsilon(1).*kz(kx,epsilon(2)));

	% if the length is two, return the two layer Fresnel, otherwise
	% recursively go through the layers
	if length(epsilon) == 2
		out = r12(kx,epsilon);
	else
		out = (r12(kx,epsilon)+nlayerfresnel(k0,kx,epsilon(2:end),d(2:end)).*exp(2.0i.*kz(kx,epsilon(2)).*d(2)))./...
					(1+r12(kx,epsilon).*nlayerfresnel(k0,kx,epsilon(2:end),d(2:end)).*exp(2.0i.*kz(kx,epsilon(2)).*d(2)));
    end
    
end
