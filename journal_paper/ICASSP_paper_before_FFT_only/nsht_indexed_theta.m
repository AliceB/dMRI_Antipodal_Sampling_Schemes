function thetas  = nsht_indexed_theta(L)
% nsht_indexed_theta - Determine sample positions along co-latutude (theta)
%
% Determine L sample positions along co-latutude (theta) for placement of
% L iso-latitude rings in the sampling scheme. The rings are placed such that the ring with more
% number of samples is placed nearer to the equator.
%
% Default usage is given by
%
%   thetas  = nsht_indexed_theta(L)
%
% where L is the harmonic band-limit , theta denotes the sample positions
% for L rings
%
%       
%
% Author: Zubair Khalid
%
% NSHT package to compute spherical harmonic transform of band-limited
% signal
% Copyright (C) 2014  Zubair Khalid
% See LICENSE.txt for license details

%%
TT_temp = pi*(2*(0:1:(L-1))+1)/(2*L-1);
[ ~, T_index] = sort(abs(TT_temp-pi/2),'descend');
thetas = TT_temp((T_index));
thetas(1) = pi - 1e-15; %so don't get theta at 0
%add to thetas the antipodal thetas
thetas = [thetas, (pi - thetas)]; 

end

