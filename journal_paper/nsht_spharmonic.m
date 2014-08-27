function [Ylm THETA PHI] = nsht_spharmonic(L,ell,m)
%% THIS FUNCTION REQUIRES AN UPDATE FOLLOWING THE nsht_legmat 
%%
% nsht_spharmonic - computes the spherical harmonic of degree L and order m
% over L^2 number of samples with the condition that
% ell < L
% Please see the paper for details about the sampling theorem



if ( ell >= L) 
  error('Require ell < L.');
end


if ( m > ell) 
  error('Require m <= ell.');
end

%[THETA, PHI] = sshtopt_sampling_points(L);

Ylm = exp(1i*m*PHI);

m=abs(m); % if m is neagtive

P = zeros(size(THETA));

Sc = zeros(size(P)); % scaling matrix contains exponents of 10


Km  = norm_factor_first(m); % find scaling factor
    
factor_remain = Km/10^floor(log10(Km));


for ii=1:1:length(THETA)
    TT = THETA(ii);

    [factor_remain2 power_10]  = sine_factor(TT, m);
    temp = floor(log10(Km)) + power_10;
    dlm = factor_remain*factor_remain2*(-1)^m*sqrt((2*m+1)/2)/2^m/sqrt(2*pi);


%%%%%%%%%%%%%%%%%%%%%


    dlm_1 = zeros(size(TT));

    for ell=m:1:ell

        if abs(dlm)<10^(-2)
            temp = temp-2;
            dlm = dlm*10^2;
            dlm_1 = dlm_1*10^2;
        end
        if abs(dlm)>10^(2)
            temp = temp+2;
            dlm = dlm/10^2;
            dlm_1 = dlm_1/10^2;
        end
        P(ii) = dlm;

        Sc(ii) =  temp;


        [ dlm1 ] = Kostelec_recursion_scaled(dlm, dlm_1,ell,m,TT );


        dlm_1 = dlm;
        dlm = dlm1;

    end


    
end

%% Now create the spharmonic and take into account phi contribution
for ii=0:length(THETA)-1;    
        Ylm(ii^2+1:(ii+1)^2) = 10.^Sc(ii+1).*  P(ii+1) *Ylm(ii^2+1:(ii+1)^2);
end


end










function [ dlm1 ] = Kostelec_recursion_scaled(dlm, dlm_1,ell,m,theta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% dlm1 dlm at ell+1
% dlm_1  ell-1
% dlm - el 
%% must call for each theta
%%

if (ell==0)
dlm1 = sqrt((2*ell+3)/(2*ell+1))* ((ell+1)*(2*ell+1)/sqrt( ( (ell+1)^2-m^2 )*((ell+1)^2)))*(cos(theta).*dlm);
else
dlm1 = sqrt((2*ell+3)/(2*ell+1))* ((ell+1)*(2*ell+1)/sqrt( ( (ell+1)^2-m^2 )*((ell+1)^2)))*(cos(theta).*dlm) - ...
       sqrt((2*ell+3)/(2*ell-1))* (sqrt((ell^2-m^2)*ell^2)/sqrt( ( (ell+1)^2-m^2 )*((ell+1)^2)))*(ell+1)*dlm_1/ell; 
end
end




function Km  = norm_factor_first(m)

Km=1;


for ii=1:1:m
    
   Km = Km*sqrt((2*m-ii+1)/(m-ii+1)); 
       
end


end


function [factor_remain power_10]  = sine_factor(TT, m)

factor_remain = 1;
power_10=0;
    for ii=1:1:m
        factor_remain = factor_remain*sin(TT);
        if abs(factor_remain)<10^(-2)
            power_10 = power_10-2;
            factor_remain = factor_remain*10^2;
        end
    end

end