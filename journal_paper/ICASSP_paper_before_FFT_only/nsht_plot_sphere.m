function [thetas phis X Y Z] = nsht_plot_sphere(f, L)
% nsht_plot_sphere -  Plots the function on the sphere over the
% sample points in the proposed sampling scheme.
%
% Default usage is
%         [thetas phis X Y Z] = nsht_plot_sphere(f, L), 
%
% L is the harmonic band-limit. f is a vector of size L^2 and a signal on the sphere
% evauated over L^2 sampling points.
%       
% Note: Since the sampling scheme does not form regular grid. This
% functions uses the matlab function TriScatteredInterp for interpolation
% of the function over the regular grid.  Due to interpolation, visual
% artifacts may be observed.
%
%
%
% Author: Zubair Khalid
%
% NSHT package to compute spherical harmonic transform of band-limited
% signal
% Copyright (C) 2014  Zubair Khalid
% See LICENSE.txt for license details
%
%%
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

[thetas, phis] = nsht_sampling_points(L);


thetas_extended = zeros(size(phis));

X = zeros(size(phis));
Y = zeros(size(phis));
Z = zeros(size(phis));

for ii=0:1:length(thetas)-1
thetas_extended(ii^2+1:(ii+1)^2) = thetas(ii+1);
end



TT = 0:pi/L:pi;
FF = 0:pi/L:2*pi;
[TT,FF] = meshgrid(TT,FF);            

F = TriScatteredInterp(thetas_extended', phis', transpose(f), 'natural');
V  = F(TT,FF);




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


