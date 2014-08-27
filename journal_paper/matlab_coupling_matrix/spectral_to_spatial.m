function [ f_spatial ] = spectral_to_spatial( f_spectral, L, L_res )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% L band-limit of the f_spectral
% L_res : spatial resolution.. number of points along theta

f_spatial=0; % initialize

for ell=0:1:L
    for m=0:1:ell
        [Ylm]=spharmonics_MW(ell,m,L_res,0);
        f_spatial = f_spatial + f_spectral(ell^2+ell+m+1)*Ylm;
        f_spatial = f_spatial + f_spectral(ell^2+ell-m+1)*(-1)^m*conj(Ylm); % -ve m
    end
end
        
            
end

