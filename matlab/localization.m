clc;
hold off;
f_location_x = fopen('Scatterers_x.txt');
location_x = fscanf(f_location_x, '%f');
f_location_y = fopen('Scatterers_y.txt');
location_y = fscanf(f_location_y, '%f');
size(location_x)
f_counting_64_00 = fopen('visting_counter_00.txt');
counting_64_00 = fscanf(f_counting_64_00, '%d');

f_counting_64_0 = fopen('visting_counter_0.txt');
counting_64_0 = fscanf(f_counting_64_0, '%d');

f_counting_64 = fopen('visting_counter_1.txt');
counting_64 = fscanf(f_counting_64, '%d');

f_counting_64_2 = fopen('visting_counter_2.txt');
counting_64_2 = fscanf(f_counting_64_2, '%d');

location_x = location_x(1:64);
location_y = location_y(1:64);

min_64_00 = min(counting_64_00);
s = counting_64_00./min_64_00*50;
c = sqrt(counting_64_00);
subplot(2,2,1);
scatter(location_x,location_y,s,c,'fill');
title('NSA = 500');

min_64_0 = min(counting_64_0);
s = counting_64_0./min_64_0*50;
c = sqrt(counting_64_0);
subplot(2,2,2);
scatter(location_x,location_y,s,c,'fill');
title('NSA = 50000');

min_64 = min(counting_64);
s = counting_64./min_64*50;
c = sqrt(counting_64);
subplot(2,2,3);
scatter(location_x,location_y,s,c,'fill');
title('NSA = 50000*10');

min_64_2 = min(counting_64_2);
s = counting_64_2./min_64_2*50;
c = sqrt(counting_64_2);
subplot(2,2,4);
scatter(location_x,location_y,s,c,'fill');
title('NSA = 50000*50');


