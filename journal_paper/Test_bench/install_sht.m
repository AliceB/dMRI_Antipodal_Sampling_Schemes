%% Define Paths
addpath ../

%% Create Mex Files for Computation of Legendre Matrix
warning off;
clc;
L=512;
Lmax = 4096*50;
m=0;
%thetas = nsht_indexed_theta(L);
%
clear mex
%codegen -args {thetas L m} nsht_legmat.m
codegen  nsht_legmat.m -args {coder.typeof(double(0),[1 Lmax],[0 1]) L m}
