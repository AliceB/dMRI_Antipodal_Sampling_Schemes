function [flm T_mat_inv] = nsht_forward(f,L)
% nsht_forward - Computes forward spherical harmonic transform
%
% Computes forward spherical harmonic transform based on the sampling scheme presented in paper
%
% Default usage is given by
%
%   f = nsht_forward(flm, L)
%
% where L is the harmonic band-limit, f is the vector of L^2 evaluated over the sampling scheme 
% and flm is the vector of L^2 harmonic coefficients
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

if ~(sum(size(f))== length(f)+1) 
      error('f must be a vector and not a matrix');
end



    if isreal(f)
        flm = zeros(1,(L)^2);
    
        [THETA, FI] = nsht_sampling_points(L); 
        
        T_mat_inv=0;


        fft_v = zeros(1,(L)^2);

        for m=L-1:-2:0
            [P_1 Sc_1] = nsht_legmat_mex(THETA, L, m);
            P_mat_1 = P_1.*10.^Sc_1;
            
            %for m-1
            if m~=0
                [P_2 Sc_2] = nsht_legmat_mex(THETA, L, m-1);
                P_mat_2 = P_2.*10.^Sc_2;
            end
            
            %uniform sampling - FFT version
            fft_v_bit = fft(f((m^2-m+2)/2:(m+1)*(m+2)/2))/length(f((m^2-m+2)/2:(m+1)*(m+2)/2));       

            fft_v(m^2+1:(m+1)^2) = fft_v_bit;

            tmatinv = tic;
            flm_t_1 = nsht_inversion( P_mat_1, fft_v,m,L);
            if m~=0
                flm_t_2 = nsht_inversion( P_mat_2, fft_v,m-1,L); %coefficients for m-1
            end
            T_mat_inv = T_mat_inv + toc(tmatinv);

            % perform spatial elimitation
            [ f ] = nsht_spatial_elimination(P_mat_1, f(1:m*(m-1)/2),flm_t_1, FI(1:m*(m-1)/2), m );
             if m~=0
                [ f ] = nsht_spatial_elimination(P_mat_2, f(1:m*(m-1)/2),flm_t_2, FI(1:m*(m-1)/2), m-1 );
             end
            
            %store coefficients calculated
            for el=m:1:L-1
               flm(el^2+el+m+1) = flm_t_1(el-m+1);
               flm(el^2+el-m+1) = (-1)^m*conj(flm_t_1(el-m+1));
            end
            
           if m~=0 
               n = m-1;
               for el=n:1:L-1

                   flm(el^2+el+n+1) = flm_t_2(el-n+1);
                   flm(el^2+el-n+1) = (-1)^n*conj(flm_t_2(el-n+1));
               end
           end
            
        end
    else % if f is complex
        flm_real = zeros(1,(L)^2);
        flm_imag = zeros(1,(L)^2);
        f_real = real(f);
        f_imag = imag(f);
        fft_v_real = zeros(1,(L)^2);
        fft_v_imag = zeros(1,(L)^2);

        [THETA, FI] = nsht_sampling_points(L); 
               
        T_mat_inv=0;

        %as zero for m odd
        for m=L-1:-2:0
            
            [P_1 Sc_1] = nsht_legmat_mex(THETA, L, m);
            P_mat_1 = P_1.*10.^Sc_1;
            %for m-1
             if m~=0
                [P_2 Sc_2] = nsht_legmat_mex(THETA, L, m-1);
                P_mat_2 = P_2.*10.^Sc_2;
             end
                 
            %FFT            
            fft_v_bit_real = fft(f_real((m^2-m+2)/2:(m+1)*(m+2)/2))/length((m^2-m+2)/2:(m+1)*(m+2)/2);
            fft_v_bit_imag = fft(f_imag((m^2-m+2)/2:(m+1)*(m+2)/2))/length(f((m^2-m+2)/2:(m+1)*(m+2)/2));
            
                      
            fft_v_real(m^2+1:(m+1)^2) = fft_v_bit_real;
            fft_v_imag(m^2+1:(m+1)^2) = fft_v_bit_imag;

            tmatinv = tic;
            flm_t_real = nsht_inversion( P_mat_1, fft_v_real,m,L);
            flm_t_imag = nsht_inversion( P_mat_1, fft_v_imag,m,L);
            if m ~= 0
                flm_t_2_real = nsht_inversion( P_mat_2, fft_v_real,m-1,L); %coefficients for m-1
                flm_t_2_imag = nsht_inversion( P_mat_2, fft_v_imag,m-1,L); %coefficients for m-1
            end
            T_mat_inv = T_mat_inv + toc(tmatinv);

            % Following needs to be updated/optimized
            f_real = nsht_spatial_elimination(P_mat_1, f_real(1:m*(m-1)/2),flm_t_real, FI(1:m*(m-1)/2), m );
            f_imag  = nsht_spatial_elimination(P_mat_1, f_imag(1:m*(m-1)/2),flm_t_imag, FI(1:m*(m-1)/2), m );
             if m~=0
                f_real = nsht_spatial_elimination(P_mat_2, f_real(1:m*(m-1)/2),flm_t_2_real, FI(1:m*(m-1)/2), m-1 );
                f_imag  = nsht_spatial_elimination(P_mat_2, f_imag(1:m*(m-1)/2),flm_t_2_imag, FI(1:m*(m-1)/2), m-1 );
             end
             
            for el=m:1:L-1
               flm_real(el^2+el+m+1) = flm_t_real(el-m+1);
               flm_real(el^2+el-m+1) = (-1)^m*conj(flm_t_real(el-m+1));
               flm_imag(el^2+el+m+1) = flm_t_imag(el-m+1);
               flm_imag(el^2+el-m+1) = (-1)^m*conj(flm_t_imag(el-m+1));

            end
            
            if m~=0
               n = m-1;
               for el=n:1:L-1

                   flm_real(el^2+el+n+1) = flm_t_2_real(el-n+1);
                   flm_real(el^2+el-n+1) = (-1)^n*conj(flm_t_2_real(el-n+1));
                   flm_imag(el^2+el+n+1) = flm_t_2_imag(el-n+1);
                   flm_imag(el^2+el-n+1) = (-1)^n*conj(flm_t_2_imag(el-n+1));
                end
            end
        end

        flm = flm_real + 1i*flm_imag;
        
    end % end reality if/else
end


function f_t = nsht_spatial_elimination(P_mat, f_t, flm_t,  FI_t, m)
% nsht_spatial_elimination - removes the m and -ve m order coefficients from the signal
% when flm is passed, assuming fl(-m) are for real signal. Input to
% function would be spatial signal 2m+1 coefficents and order m, TT and FF. The output
% would be spatial signal with coefficients removed.

% f_t = truncated signal 
% flm_t spectrum of order m contains L-m+1 coefficients, the first one
% corresponds to degree el=m and the last one for degree el=L-1

% THETA_t, FI_t, the truncated spatial grid.
% m order
% L is the band-limit
f_temp=zeros(size(FI_t));
f_temp_neg=zeros(size(FI_t));

if mod(m,2) == 0 % m even
    gm = transpose(P_mat(:,1:m/2))*flm_t; 
    gm_neg = conj(transpose(P_mat(:,1:m/2))*flm_t);
    %assigning gm(theta) to the correct  f(theta, phi)
    for ii=0:m/2-1
       f_temp(2*ii^2-ii+1:2*ii^2+3*ii+1) = gm(ii+1); 
       f_temp_neg(2*ii^2-ii+1:2*ii^2+3*ii+1) = gm_neg(ii+1); 
    end
else
     gm = transpose(P_mat(:,1:(m+1)/2))*flm_t; 
    gm_neg = conj(transpose(P_mat(:,1:(m+1)/2))*flm_t);
    %assigning gm(theta) to the correct  f(theta, phi)
    for ii=0:(m+1)/2-1
       f_temp(2*ii^2-ii+1:2*ii^2+3*ii+1) = gm(ii+1); 
       f_temp_neg(2*ii^2-ii+1:2*ii^2+3*ii+1) = gm_neg(ii+1); 
    end
end




% no need of condition on m here because the spatial elimination is not run
% for m=0. (Verified. 15/11/2013)

f_t = f_t - f_temp.*exp(1i*m*FI_t) - f_temp_neg.*exp(-1i*m*FI_t);

end


function flm_t  = nsht_inversion(P_mat, fft_v,m,L)
% nsht_inversion - when passed with fft vector and order m,
% band-limit L, returns the spherical harmonic coefficients. This function
% would perform matrix inversion.

% fft_v should have the same structure as of FI
% flm_t is the vector which contains spectral coefficeints [fm^m, f(m+1)^m
% fL^m];

% Ax = b, first form b from fft_vec

b = zeros(1,L-m); 
for el = m:L-1
      b(el-m+1) = fft_v(m+el^2+1);     
end

%keep non-zero values
if mod(m,2) == 0 % m even
    b = b(1:2:L-m);
    P = P_mat(1:2:end,m/2+1:end);
else
    b = b(2:2:L-m);
    P = P_mat(2:2:end,(m+1)/2+1:end);
end

global conditionNum;
conditionNum(m + 1)= cond(P);

flm_short = (transpose(P))\transpose(b);
flm_t = zeros(1,L-m);
if mod(m,2) == 0 % m even
    for i = 1: length(flm_short)
        flm_t((i*2)-1) = flm_short(i);
    end
else
     for i = 1: length(flm_short)
        flm_t((i*2)) = flm_short(i);
    end

end

flm_t = transpose(flm_t);

end

