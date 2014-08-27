L = 40; % band-limit. I think the maximum band-limit supported is L=85 because of matlab limit of factorial

theta_c = pi/4; % central angle of the localized spatial region 

m = 0; % order m = 0,1,2..L

% The Gll_mat is the coupling matrix for m-order spherical harmonic
% coefficients. Size: (L-m+1)x(L-m+1)

% These are computed exactly using the analytic formulation and using Wigner 3j
% symbols which are alculated using the routine by Simons and are verified
% through direct computation.

[ Gll_mat ] = Gll_matrix_order(L+1, m,theta_c ); 

%% Eigen-decomposition

[V,D] = eig(Gll_mat);
eig_v = flipud(real(diag(D)));
V  = fliplr(V); % just flipping so that most concentarted is at index 1

figure;
plot(eig_v);

%% Computation of eigenfunctions in the spatial domain
ii=3; % ii-th eigenfunction. ii = 0,1...(L-m+1)

L_res = 128; % choosing spatial resolution

f_spectral = zeros(1,(L+1)^2);

for ell=m:1:L
    f_spectral(ell^2+ell+m+1) = V(ell-m+1,ii);
end

% Note that the coupling matrix is computed only for positive m. For negative m, the eigenfunction would be just a complex conugate
% of eigenfucntion for positive m.

[ f_spatial ] = spectral_to_spatial( f_spectral, L, L_res );

%% Plotting
[ THETA,FI, X,Y,Z ] = slsht_sampling_points(L_res);

eig_value_for_eigenfunction = eig_v(ii)
figure('Color', [1 1 1]);
h = surf(X,Y,Z,abs(f_spatial)); % can be changed to imaginary
set(h, 'LineStyle', 'none');
light; lighting phong;
view([-1,-1,1]);
axis equal;
axis([-1.001 1.001 -1.001 1.001 -1.001 1.001])
h_bar = colorbar('location','EastOutside');
