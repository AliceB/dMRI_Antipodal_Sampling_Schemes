function [THETA, FI] = nsht_sampling_points(L)
% nsht_sampling_points - Determine the sampling points on the sphere for
% the sampling scheme proposed in the following paper.
% 
% Paper
% 
% Default usage is given by
%
%   [THETA, FI] = nsht_sampling_points(L)
%
% where L is the harmonic band-limit. THTEA denotes the \theta vector of
% size L. PHI contains the position of samples along longitude and is a
% vector of size L^2.
%       
% Author: Zubair Khalid
%
% NSHT package to perform spherical harmonic transforms

% Check arguments.
% Check arguments.


THETA = nsht_ordered_theta(L);

FI = zeros(1,L*(L+1)/2);

    %%
    start=0;
    for LL=0:2:L-1
       FF_temp = pi*(2*(0:1:2*LL))/(2*LL+1); % 
       FI(start+1:start+length(FF_temp)) = FF_temp;
      start = start+length(FF_temp); 

    end

  
end



