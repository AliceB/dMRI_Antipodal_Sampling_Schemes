function flm = nsht_forward(f,L)
% nsht_forward - Computes forward spherical harmonic transform
%
% Computes forward spherical harmonic transform based on the sampling
% scheme presented in paper.
%
% Default usage is given by
%
%   flm = nsht_forward(flm, L)
%
% where L is the harmonic band-limit, f is the vector of L^2 evaluated over the sampling scheme 
% and flm is the vector of L^2 harmonic coefficients
%
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

% Check arguments.
if ~isreal(L)
      error('Harmonic band-limit must be real');
end

if ~(sum(size(f))== length(f)+1) 
      error('f must be a vector and not a matrix');
end

if ~(length(f)==(L+1)*L/2)
      error('f must be a vector of size (L+1)*L/2');
end

    if isreal(f)
        flm = zeros(1,(L)^2);

        [THETA, FI] = nsht_sampling_points(L); % make sure this returns the re-ordered THETA



        fft_v = zeros((L)^2);

        for m=L-1:-1:0
            [P Sc] = nsht_legmat_mex(THETA, L, m);
            P_mat = P.*10.^Sc;
                    
            
            %use NDFT as phi doesn't start at zero for m odd
            if mod(m,2)==0
                phi = FI(0.5*m^2-0.5*m+1:0.5*m^2+1.5*m+1).';
                fbit = f(0.5*m^2-0.5*m+1:0.5*m^2+1.5*m+1).';
                E = exp(1i*(phi)*(-m:1:m));
            else 
                n = m +1;
                phi= mod(FI(0.5*n^2-0.5*n+1:0.5*n^2+1.5*n+1)+ pi, 2*pi).';  %antipodal symmetry
                fbit = f(0.5*n^2-0.5*n+1:0.5*n^2+1.5*n+1).';
                E = exp(1i*(phi)*(-n:1:n));
            end
            
            
            fft_v_bit = (inv(E)*fbit).';    

             
            if mod(m,2)==0
                fft_v_bit = ifftshift(fft_v_bit);
                fft_v(m^2+1:(m+1)^2) = fft_v_bit;
     
            else
                fft_v_bit = ifftshift(fft_v_bit(2:(end-1)));
                fft_v(m^2+1:(m+1)^2) = fft_v_bit;
                
            end


            flm_t = nsht_inversion( P_mat, fft_v,m,L);
           
            if mod(m,2)==0
                f = nsht_spatial_elimination(P_mat, f(1:0.5*(m^2+3*m+2)),flm_t, FI(1:0.5*(m^2+3*m+2)), m );
            else
                [ f ] = nsht_spatial_elimination(P_mat, f(1:0.5*m*(m+1)),flm_t, FI(1:0.5*m*(m+1)),  m);               
            end

            % Following needs to be updated/optimized
            
            for el=m:1:L-1

               flm(el^2+el+m+1) = flm_t(el-m+1);
               flm(el^2+el-m+1) = (-1)^m*conj(flm_t(el-m+1));
            end
        end
    else % if f is complex
        flm_real = zeros(1,(L)^2);
        flm_imag = zeros(1,(L)^2);
        f_real = real(f);
        f_imag = imag(f);
        fft_v_real = zeros(1,(L)^2);
        fft_v_imag = zeros(1,(L)^2);

        [THETA, FI] = nsht_sampling_points(L); % make sure this returns the re-ordered THETA




        for m=L-1:-1:0
            [P Sc] = nsht_legmat_mex(THETA, L, m);
            P_mat = P.*10.^Sc;

            
             %use NDFT as phi doesn't start at zero for m odd
            if mod(m,2)==0
                phi = FI(0.5*m^2-0.5*m+1:0.5*m^2+1.5*m+1).';
                fbitReal = f_real(0.5*m^2-0.5*m+1:0.5*m^2+1.5*m+1).';
                fbitImag = f_imag(0.5*m^2-0.5*m+1:0.5*m^2+1.5*m+1).';
                E = exp(1i*(phi)*(-m:1:m));
            else 
                n = m +1;
                phi= mod(FI(0.5*n^2-0.5*n+1:0.5*n^2+1.5*n+1)+ pi, 2*pi).';  %antipodal symmetry
                fbitReal = f_real(0.5*n^2-0.5*n+1:0.5*n^2+1.5*n+1).';
                fbitImag = f_imag(0.5*n^2-0.5*n+1:0.5*n^2+1.5*n+1).';
                E = exp(1i*(phi)*(-n:1:n));
            end
            
            
            fft_v_bit_real = (inv(E)*fbitReal).';            
            fft_v_bit_imag = (inv((E))*fbitImag).';
           
            
            if mod(m,2)==0
                fft_v_bit_real = ifftshift(fft_v_bit_real);
                fft_v_bit_imag = ifftshift(fft_v_bit_imag);
                fft_v_real(m^2+1:(m+1)^2) = fft_v_bit_real;
                fft_v_imag(m^2+1:(m+1)^2) = fft_v_bit_imag;
            else
                fft_v_bit_real = ifftshift(fft_v_bit_real(2:(end-1)));
                fft_v_bit_imag = ifftshift(fft_v_bit_imag(2:(end-1)));
                fft_v_real(m^2+1:(m+1)^2) = fft_v_bit_real;
                fft_v_imag(m^2+1:(m+1)^2) = fft_v_bit_imag;
                
            end
           


            flm_t_real = nsht_inversion( P_mat, fft_v_real,m,L);
            flm_t_imag = nsht_inversion( P_mat, fft_v_imag,m,L);

            if mod(m,2)==0
                f_real = nsht_spatial_elimination(P_mat, f_real(1:0.5*(m^2+3*m+2)),flm_t_real, FI(1:0.5*(m^2+3*m+2)), m );
                f_imag  = nsht_spatial_elimination(P_mat, f_imag(1:0.5*(m^2+3*m+2)),flm_t_imag, FI(1:0.5*(m^2+3*m+2)), m );
            else
                f_real = nsht_spatial_elimination(P_mat, f_real(1:0.5*m*(m+1)),flm_t_real, FI(1:0.5*m*(m+1)), m );
                f_imag  = nsht_spatial_elimination(P_mat, f_imag(1:0.5*m*(m+1)),flm_t_imag, FI(1:0.5*m*(m+1)), m );
            end
            
            for el=m:1:L-1
               flm_real(el^2+el+m+1) = flm_t_real(el-m+1);
               flm_real(el^2+el-m+1) = (-1)^m*conj(flm_t_real(el-m+1));
               flm_imag(el^2+el+m+1) = flm_t_imag(el-m+1);
               flm_imag(el^2+el-m+1) = (-1)^m*conj(flm_t_imag(el-m+1));

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

if mod(m,2)==0
     gm = (transpose(P_mat(:,1:(m+1)))*flm_t);
    gm_neg = conj(transpose(P_mat(:,1:(m+1)))*flm_t);
else
    gm = (transpose(P_mat(:,1:m))*flm_t);
    gm_neg = conj(transpose(P_mat(:,1:m))*flm_t);

end

f_temp=zeros(size(FI_t));
f_temp_neg=zeros(size(FI_t));

    for ii=0:2:m
        f_temp(0.5*ii^2-0.5*ii+1:0.5*ii^2+1.5*ii+1) = gm(ii+1); 
        f_temp_neg(0.5*ii^2-0.5*ii+1:0.5*ii^2+1.5*ii+1) = gm_neg(ii+1); 
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
for el=m:1:L-1
      b(el-m+1) = fft_v(m+el^2+1);     
end


P = P_mat(:,m+1:end);


flm_t = (transpose(P))\transpose(b);



end

