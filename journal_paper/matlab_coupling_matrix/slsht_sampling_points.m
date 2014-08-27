function [ THETA,FI, X,Y,Z ] = slsht_sampling_points(L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


alpha = 0:2*pi/(2*L+1):2*pi- 2*pi/(2*L+1);
beta = 0:2*pi/(2*L):2*pi;

beta = beta(1:L+1);

[FI,THETA] = meshgrid(alpha,beta);

[X, Y, Z] = slsht_s2c(THETA, FI);

end

function [X, Y, Z] = slsht_s2c(THETA, FI)

X = sin(THETA).*cos(FI);
Y = sin(THETA).*sin(FI);
Z = cos(THETA);

end
