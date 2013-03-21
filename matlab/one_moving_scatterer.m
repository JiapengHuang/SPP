clc;
clear;
close all;
hold off;
NSCAT = 65;
VAR = 100;
location_x = load('Scatterers_x_moving.txt');
location_y = load('Scatterers_y_moving.txt');

location_x = reshape(location_x,NSCAT,VAR);
size(location_x);
location_y = reshape(location_y,NSCAT,VAR);
size(location_y);

counting_65 = load('visting_counter_moving.txt');

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
end
%movie(A,10,3);
movie2avi(A,'localization.avi','compression','None');
