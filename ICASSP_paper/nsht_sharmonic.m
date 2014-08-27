function [Ylm thetas phis] = nsht_sharmonic(L,ell,m)
% nsht_sharmonic - Compute spherical harmonic over optimal sampling spatial grid
%
% Computes spherical harmonic of degree ell and order m
% over the optimal sampling grid of L^2 samples defined by L. Please see the paper for
% details about the sampling scheme. 
%
% Usage is given by
%
%  Ylm = nsht_sharmonic(ell,m,thetas,phis)
%
% where ell is degree and m is order and thetas and phis 
% are vectors of same size. Ylm is evaluated over (thetas,phis) pairs
%
% Notes:
%  - Kostelec recusrsion is implemented to compute the scaled legendre
%  coefficients for each theta \in thetas.
%  - It uses nsht_sharmonic_g which computes spherical harmonic of degree
%  over given sample points.
%
%
% Author: Zubair Khalid
%
% NSHT package to compute spherical harmonic transform of band-limited
% signal
% Copyright (C) 2014  Zubair Khalid
% See LICENSE.txt for license details
%%

% Check arguments. 
if ( ell >= L) 
  error('Require ell < L.');
end
if ( m > ell) 
  error('Require m <= ell.');
end

[thetas, phis] = nsht_sampling_points(L);

% Extend thetas to match the dimension of phis
thetas_ext_L = zeros(size(phis));
for ii=0:1:length(thetas)-1
   thetas_ext_L((ii)^2+1: (ii+1)^2) = thetas(ii+1);
end

% call nsht_sharmonic_g which computes spherical harmonic of degree ell and order m
% over given thetas_ext_L and phis points

Ylm = nsht_sharmonic_g(ell,m,thetas_ext_L,phis);


end