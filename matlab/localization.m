clc;
hold off;
location_x     = load('Scatterers_x.txt');
location_y     = load('Scatterers_y.txt');
size(location_x)
counting_64_00 = load('visting_counter_00.txt');
counting_64_0  = load('visting_counter_0.txt');
counting_64    = load('visting_counter_1.txt');
counting_64_2  = load('visting_counter_2.txt');

location_x = location_x(1:64);
location_y = location_y(1:64);

min_64_00 = min(counting_64_00);
s = counting_64_00./min_64_00*50;
c = sqrt(counting_64_00);
subplot(2,2,1);
% original matlab code, changed so it displays correctly in octave
% scatter(location_x,location_y,s,c,'fill');
plot(location_x,location_y,'o');
title('NSA = 500');

min_64_0 = min(counting_64_0);
s = counting_64_0./min_64_0*50;
c = sqrt(counting_64_0);
subplot(2,2,2);
plot(location_x,location_y,'o');
title('NSA = 50000');

min_64 = min(counting_64);
s = counting_64./min_64*50;
c = sqrt(counting_64);
subplot(2,2,3);
plot(location_x,location_y,'o');
title('NSA = 50000*10');

min_64_2 = min(counting_64_2);
s = counting_64_2./min_64_2*50;
c = sqrt(counting_64_2);
subplot(2,2,4);
plot(location_x,location_y,'o');
title('NSA = 50000*50');


