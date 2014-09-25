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
phis = zeros(1,L^2);

%%
start=0;
for LL=0:1:L-1
    FF_equispaced = pi*(2*(0:1:2*LL))/(2*LL+1); %where samples would be for the optimal dimesionality sampling scheme
    
%     %case L even
%     if mod(L,2) == 0
%         if mod(LL,2) == 0    
%             FF_available = mod(pi*(2*(0:1:2*(LL+1)))/(2*(LL+1)+1) + pi, 2*pi); %because of antipodal nature of samples      
%             FF_temp = zeros(1, 2*LL + 1);
%             for i=1:length(FF_equispaced)
%                [~,INDEX] = min(abs(FF_available - FF_equispaced(i)));   %find the closest phi sample 
%                FF_temp(i) =FF_available(INDEX);
%             end
% 
%             %%check no repeated points
%            if ~(length(unique(FF_temp)) == length(FF_temp))
%                disp('PROBLEM - PHI VALUES REPEATED');
%            end
% 
%         else
%             FF_temp = FF_equispaced;
%         end
%     else 
%case L odd
        if mod(LL,2) ~= 0    
            n = LL + 1;
            FF_available = mod(pi*(2*(0:1:2*n))/(2*n+1) + pi, 2*pi); %because of antipodal nature of samples      
            FF_temp = zeros(1, 2*LL + 1);
            for i=1:length(FF_equispaced)
               [~,INDEX] = min(abs(FF_available - FF_equispaced(i)));   %find the closest phi sample 
               FF_temp(i) =FF_available(INDEX);
            end

            %%check no repeated points
           if ~(length(unique(FF_temp)) == length(FF_temp))
               disp('PROBLEM - PHI VALUES REPEATED');
           end

        else
            FF_temp = FF_equispaced;
        end
    
%end
   phis(start+1:start+length(FF_temp)) = FF_temp;
  start = start+length(FF_temp);     
end


end
