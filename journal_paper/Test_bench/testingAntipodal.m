
install_sht;

% Generate random signal
 L = 64;
 ft = zeros(1,L^2);
 index = [1,0,5,0,17,0,37];
 for l = 1:2:7
     ft(index(l):(index(l) +(2*l-2) )) =  rand(1,(2*l-1)) + 1i*rand(1,(2*l-1));
 end
 

%% inverse transform 
fr  = nsht_inverse(ft,L);

%%
nsht_plot_sphere(real(fr), L)
nsht_plot_sphere(imag(fr), L)