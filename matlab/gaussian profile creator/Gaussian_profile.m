clc;
close all;
clear;
k_spp =  1.0541e+07;
k = k_spp + (0.016*k_spp).*randn(1,50000);
fid = fopen('Gaussian_k.txt', 'wt');
fprintf(fid, '%f\n', k);
fclose(fid);
hist(k,100);