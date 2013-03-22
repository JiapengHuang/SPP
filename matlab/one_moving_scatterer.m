clc;
clear;
close all;
hold off;
NSCAT = 65;
VAR = 100;
f_location_x = fopen('Scatterers_x_moving.txt');
location_x = fscanf(f_location_x, '%f');
f_location_y = fopen('Scatterers_y_moving.txt');
location_y = fscanf(f_location_y, '%f');

location_x = reshape(location_x,NSCAT,VAR);
size(location_x);
location_y = reshape(location_y,NSCAT,VAR);
size(location_y);

f_counting_65 = fopen('visting_counter_moving.txt');
counting_65 = fscanf(f_counting_65,'%d');

counting_65 = reshape(counting_65,NSCAT,VAR);
size(counting_65);

A = moviein(VAR);

counting = counting_65(:,1);
min_65 = min(counting);
for j = 1:VAR
    counting = counting_65(:,j);
    min_65 = min(counting);
    s = counting./min_65*300;
    c = sqrt(counting./min_65);
    figure;
    scatter(location_x(:,j),location_y(:,j), s,c,'fill');
    title('**');
    A(:,j) = getframe;
    im = frame2im(A(:,j));
    [imind,cm] = rgb2ind(im,256);
    if j == 1;
    imwrite(imind,cm,'localizations','gif', 'Loopcount',inf);
    else
    imwrite(imind,cm,'localizations','gif','WriteMode','append');
    end
end
%movie(A,10,3);
movie2avi(A,'localization.avi','compression','None');