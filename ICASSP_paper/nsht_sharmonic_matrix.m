function Ylm_mat = nsht_sharmonic_matrix(L,thetas,phis)
% nsht_sharmonic_matrix - Computes spherical harmonic of all degrees and 
% orders less than L over given spatial grid.
%
% Usage is given by
%
%  Ylm_mat = nsht_sharmonic_matrix(L,thetas,phis)
%
% where thetas and phis are vectors of same size. 
%
% Important Note: The size of the matrix is (L^2 x length(thetas))
% If the samples (thetas,phis) are taken from the optimal sampling scheme,
% the maximum size available to Matlab environment limits the maximum
% band-limit L, up to which this function can be executed. This requirement is explicitely
% checked below. 
%
% Notes:
%  - Kostelec recusrsion is implemented to compute the scaled legendre
%  coefficients for each theta \in thetas

% Check arguments. 
if ~(length(thetas)==length(phis))
    error('thetas and phis must be of same size');
end

if ~(sum(size(thetas))== length(thetas)+1) 
      error('thetas must be a vector and not a matrix');
end

if ~(sum(size(phis))== length(phis)+1) 
      error('phis must be a vector and not a matrix');
end

Mem_status = memory;
Mem_max = Mem_status.MemAvailableAllArrays;
if (length(thetas)*L^2)>= Mem_max
          error('phis must be a vector and not a matrix');
end
%%

tic;
Ylm_mat = zeros(L^2,length(thetas));
toc;

%%

% find leg_mat over unique thetas in thetas
thetas_unique = unique(thetas);

%%
for m=0:1:L-1
    [P Sc] = nsht_legmat_mex(thetas_unique, L, abs(m));
    P_mat = P.*10.^Sc;
    
    for ell=m:1:L-1
        for ii=1:length(thetas)
            Ylm_mat(ell^2+ell+m+1,ii) = P_mat(ell-abs(m)+1,(thetas_unique==thetas(ii)));
        end
    end

    for ell=m:1:L-1
        Ylm_mat(ell^2+ell+m+1,:) = Ylm_mat(ell^2+ell+m+1,:).*exp(1i*m*phis);
    end
end
%%
% for negative orders

for ell=1:1:L-1
    for m=1:1:ell
        Ylm_mat(ell^2+ell-m+1,:)  = (-1)^m*conj(Ylm_mat(ell^2+ell+m+1,:));    
    end
end
%%
end

