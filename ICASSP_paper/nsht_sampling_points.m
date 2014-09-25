function [thetas, phis] = nsht_sampling_points(L)
% nsht_sampling_points - Determine the sampling points on the sphere for
% the sampling scheme proposed in the following paper.
% 
% Paper
% 
% Default usage is given by
%
%   [thetas, phis] = nsht_sampling_points(L)
%
% where L is the harmonic band-limit. thetas denotes the optimal \theta vector of
% size L. phis contains the position of samples along longitude and is a
% vector of size L^2. See the paper for details.
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
thetas = nsht_ordered_theta(L);
phis = [];

%%
start=0;
for LL=0:2:L-1       
  FF_temp = pi*(2*(0:1:2*LL))/(2*LL+1);
  phis = [phis, FF_temp];
  start = start+length(FF_temp);     
end


end
