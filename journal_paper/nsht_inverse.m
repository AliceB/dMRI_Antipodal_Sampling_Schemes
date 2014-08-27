function [f ] = nsht_inverse(flm,L)
% nsht_inverse - Computes inverse spherical harmonic transform
%
% Computes inverse spherical harmonic transform based on the sampling scheme presented in paper
%
% Default usage is given by
%
%   f = nsht_inverse(flm, L, <options>)
%
% where L is the harmonic band-limit, flm is the vector of L^2
% harmonic coefficients and f is the vector of L^2 evaluated over the sampling scheme
%
%
% Option consist of parameter type and value pairs.  Valid options
% include:
%  'Reality'         = { false        [do not assume f real (default)],
%       
% Author: Zubair Khalid
%
% NSHT package to perform spherical harmonic transforms

% Check arguments.
if ~isreal(L)
      error('Harmonic band-limit must be real');
end

if ~(sum(size(flm))== length(flm)+1) 
      error('flm must be a vector and not a matrix');
end


if ~(length(flm)==L^2)
      error('flm must be a vector of size L^2');
end


[THETA, FI] = nsht_sampling_points(L);
    
    % initialize
f = zeros(size(FI));

% iteratvely compute the contribution for different orders m and -m
for m=0:1:L-1    
    
    [P Sc] = nsht_legmat_mex(THETA, L, m);
    fm = zeros(1,L-m);
    fm_neg = zeros(1,L-m);


    for el=m:1:L-1
        fm(el-m+1) = flm(el^2+el+m+1);
        fm_neg(el-m+1) = flm(el^2+el-m+1);
    end
    
    P_mat = P.*10.^Sc;
    gm = fm*P_mat;
    gm_neg = (-1)^m*fm_neg*P_mat;

    f_temp=zeros(size(FI));
    f_temp_neg=zeros(size(FI));
    
    
    for ii=0:length(THETA)-1
       f_temp(2*ii^2-ii+1:2*ii^2+3*ii+1) = gm(ii+1); 
       f_temp_neg(2*ii^2-ii+1:2*ii^2+3*ii+1) = gm_neg(ii+1); 
    end

    
    if m==0
        f = f+ f_temp;
    else
         f = f+ f_temp.*exp(1i*m*FI) + f_temp_neg.*exp(-1i*m*FI);
    end

end
        
end

