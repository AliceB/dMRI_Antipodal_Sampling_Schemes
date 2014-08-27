function [ w ] = reflected_wts( mm )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% different from jason paper, here L_h is L-1 or L = L_h+1


if mod(mm,2)==0
   w = 2./(1-mm^2) ;
else
    w=0;
end



if mm==1
   w = 1i*pi/2 ;
end

if mm==-1
   w = -1i*pi/2 ;
end



end

