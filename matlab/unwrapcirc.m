% unwrapcirc.m
% takes a cone image and, given the radius, plots the cone as a function of
% angle along with a histogram of the pixel intensities.

% filename of cone to unwrap
a = imread('out_s2000_w7_5000_50000.png');

% size of image
dims = size(a);

% approximate radius
r = 206;

% number of points to use for unwraping
N = ceil(2*pi*r);

unwrappedcirc = zeros(1,N);

% unwrap 
i=1;
for t = linspace(0,2*pi,N)
	unwrappedcirc(i) = a(ceil(r.*cos(t) + dims(1)./2), ceil(r.*sin(t) + dims(2)./2));
	i=i+1;
end

% normalize unwrappedcirc
unwrappedcirc = unwrappedcirc./max(unwrappedcirc(:));

t = linspace(0,2*pi,N);
figure(1);
polar(t,1-unwrappedcirc);
xlabel('radial angle');
ylabel('intensity');

figure(2);
hist(unwrappedcirc);

