function [X_vec Y_vec Z_vec ] = nsht_plot_sampling(L)
% nsht_plot_sampling -  Plots the sampling points on the sphere over the
% surface of the sphere. 
%
% Default usage is
%         [X_vec Y_vec Z_vec ] = nsht_plot_sampling(L), 
%
% L is the harmonic band-limit. X_vec, Y_vec, Z_vec are the vectors, each
% of length L^2 containing the x, y, z coordinates of the samples
% respectively.
%       
% Author: Zubair Khalid
%
% NSHT package to perform spherical harmonic transforms
%
% Check arguments.
if ~isreal(L)
      error('Harmonic band-limit must be real');
end

%parameters
VAR = 0; %uniform sampling along phi
LL=80;

[x y z] = sphere(LL);

figure('Color', [1 1 1]);
colormap(gray);

h = surf(0.99*x,0.99*y,0.99*z,ones(size(x)));
set(h,'Linestyle', 'none');

hold on;
axis off;
axis equal;


[THETA, FI] = nsht_sampling_points(L, VAR);


X_vec = zeros(size(FI));
Y_vec = zeros(size(FI));
Z_vec = zeros(size(FI));

%%


%% longitudes

TT = [pi/5 2*pi/5 3*pi/5 4*pi/5];
FF = 0:0.01:2*pi;

for i=0:1:length(TT)-1
        Xm = 0.99*sin(TT(i+1)).*cos(FF);
        Ym = 0.99*sin(TT(i+1)).*sin(FF);
        Zm = 0.99*cos(TT(i+1)).*ones(size(FF));
        plot3(Xm,Ym,Zm, '--','Color',[0 0 0], 'linewidth',0.5);
end
%% latitudes
FF = pi/4:pi/4:2*pi
TT = 0:0.01:pi;

for i=0:1:length(FF)-1
        Xm = 0.99*sin(TT).*cos(FF(i+1));
        Ym = 0.99*sin(TT).*sin(FF(i+1));
        Zm = 0.99*cos(TT);
        plot3(Xm,Ym,Zm, '--','Color',[0 0 0], 'linewidth',0.5);
        hold on;
end
% 

%%


temp=1;
for i=0:1:length(THETA)-1
    for j=i^2+1:1:(i+1)^2
        Xm = 1*sin(THETA(i+1)).*cos(FI(j));
        Ym = 1*sin(THETA(i+1)).*sin(FI(j));
        Zm = 1*cos(THETA(i+1));
        X_vec(temp) = Xm;
        Y_vec(temp) = Ym;
        Z_vec(temp) = Zm;
        temp=temp+1;
        plot3(Xm,Ym,Zm, '.', 'markersize', 20,  'MarkerFaceColor',[0 0 .8]);
    end
end


end

