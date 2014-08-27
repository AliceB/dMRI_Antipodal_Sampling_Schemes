function nsht_plot_spectral(flm, L)
% nsht_plot_spectral -  Plots the spectral representation of the signal as 
% a surface plot for different degress and orders.
% surface of the sphere. 
%
% Default usage is
%         nsht_plot_spectral(flm, L);
%
% flm must be real. L is the harmonic band-limit. flm is a vector of size L^2 and contains 
% the spherical harmonic coefficients of the signal with linear indexing.
%       
%
% Author: Zubair Khalid
%
% NSHT package to perform spherical harmonic transforms
%
% Check arguments.

% argument check
if ~(length(flm)==L^2)
      error('flm must be a vector of size L^2');
end

if ~isreal(flm)
      error('flm must be real');
end

x = 0:L-1;
y = -(L-1):(L-1);
clear NaN; % if variable NaN is in workspace

Flm_mat = NaN(length(x),length(y));

[x,y] = meshgrid(x,y);

for ell=0:1:L-1
    for m=-ell:1:ell
       Flm_mat(ell+1, L+m) =  flm(ell^2+ell+m+1);
    end
end


hfig = figure('Color', [1 1 1]);
set(hfig, 'Position', [200, 100, 500, 700]);
h =   surf(x,y,Flm_mat');    
set(h,'Linestyle', 'none');
shading interp
lighting phong
axis([0 1.1*L -1.1*L 1.1*L]); 
view([ 0 0 1])
set(gca,'FontSize',18)
xlabel('degree, l');
ylabel('order, m');


end
