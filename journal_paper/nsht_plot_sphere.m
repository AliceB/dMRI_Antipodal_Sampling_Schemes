function [THETA PHI X Y Z] = nsht_plot_sphere(f, L)
% nsht_plot_sphere -  Plots the function on the sphere over the
% sample points in the proposed sampling scheme.
%
% Default usage is
%         [THETA PHI X Y Z] = nsht_plot_sphere(f, L), 
%
% L is the harmonic band-limit. f is a vector of size L^2 and a signal on the sphere
% evauated over L^2 sampling points.
%       
% Note: Since the sampling scheme does not form regular grid. This
% functions uses the matlab function TriScatteredInterp for interpolation
% of the function over the regular grid.
%
%
% Author: Zubair Khalid
%
% NSHT package to perform spherical harmonic transforms
%
% Check arguments.
if ~isreal(L)
      error('Harmonic band-limit must be real');
end

if ~(sum(size(f))== length(f)+1) 
      error('f must be a vector and not a matrix');
end

if ~(length(f)==L^2)
      error('f must be a vector of size L^2');
end

if ~isreal(f)
      error('f must be real');
end

[THETA, PHI] = nsht_sampling_points(L);


THETA_extended = zeros(size(PHI));

X = zeros(size(PHI));
Y = zeros(size(PHI));
Z = zeros(size(PHI));

for ii=0:1:length(THETA)-1
THETA_extended(ii^2+1:(ii+1)^2) = THETA(ii+1);
end



TT = 0:pi/L:pi;
FF = 0:pi/L:2*pi;
[TT,FF] = meshgrid(TT,FF);            

F = TriScatteredInterp(THETA_extended', PHI', transpose(f), 'natural');
V  = F(TT,FF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%
%V = griddata(THETA_extended,PHI,f, TT,FF, 'cubic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



X = 1*sin(TT).*cos(FF);
Y = 1*sin(TT).*sin(FF);
Z = cos(TT);
        
figure('Color', [1 1 1]);
h =   surf(X,Y,Z, real(V));    
set(h,'Linestyle', 'none');
shading interp
lighting phong
axis tight
axis equal
camzoom(1.1)
set(gca,'XTick',[-1 0 1]);
set(gca,'YTick',[-1 0 1]);
set(gca,'ZTick',[-1 0 1]);
set(gca,'FontSize',18)
end


