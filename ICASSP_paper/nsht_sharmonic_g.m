function Ylm = nsht_sharmonic_g(ell,m,thetas,phis)
% nsht_sharmonic_g - Compute spherical harmonic over optimal sampling 
% spatial grid
%
% Computes spherical harmonic of degree ell and order m
% over given thetas and phis points
%
% Usage is given by
%
%  Ylm = nsht_sharmonic_g(ell,m,thetas,phis)
%
% where ell is degree and m is order and thetas and phis 
% are vectors of same size. Ylm is evaluated over (thetas,phis) pairs
%
% Notes:
%  - Kostelec recusrsion is implemented to compute the scaled legendre
%  coefficients for each theta \in thetas

%
% Author: Zubair Khalid
%
% NSHT package to compute spherical harmonic transform of band-limited
% signal
% Copyright (C) 2014  Zubair Khalid
% See LICENSE.txt for license details
%%
% Check arguments. 
if ( abs(m) > ell) 
  error('Require m <= ell.');
end

if ~(length(thetas)==length(phis))
    error('thetas and phis must be of same size');
end

Ylm = zeros(size(thetas));

% find leg_mat over unique thetas in thetas
thetas_unique = unique(thetas);
[P Sc] = nsht_legmat_mex(thetas_unique, ell+1, abs(m));
P_mat = P.*10.^Sc;

for ii=1:length(thetas)
    Ylm(ii) = P_mat(ell-abs(m)+1,find(thetas_unique==thetas(ii)));
end

if m>=0
    Ylm = Ylm.*exp(1i*m*phis);
else
    Ylm = (-1)^m*Ylm.*exp(1i*m*phis);
end

end

