%testBench.m generates reconstruction error plots for antipodal signal used
%in ICASSP paper
%author: Alice Bates
%date: September 2014

% Install SHT 
install_sht

MAX_Errors = zeros(1,13);
MEAN_Errors = zeros(1,13);


index = 1; 
for L = 1:2:25
    for j =1:10 %repeat experiment 10 times for each L


        %% Generate random antipodal signal
        flmt = zeros(1,L^2);
        for l = 1:2:L
             flmt(l^2-2*l+2:l^2) =  rand(1,(2*l-1)) + 1i*rand(1,(2*l-1));
        end
%          load('rand_L5_flm.mat')
        %% inverse transform
        ft  = nsht_inverse(flmt,L);
        
        %% forward transform
   %     tStartfwd = tic;
        flmr = nsht_forward(ft,L);
  %      toc(tStartfwd);
        
        %% Absolute error
        Error_max = max(abs(flmt-flmr));

        %%MEAN error
        Error_mean = sum(abs(flmt-flmr))/L^2;
        
        %store average RMS and MAX error for plotting
        MAX_Errors(index) = MAX_Errors(index) +  Error_max*0.1;
        MEAN_Errors(index) = MEAN_Errors(index) + Error_mean*0.1;

   end
    
    index = index + 1;
end



%% plot max and mean error
x = 1:2:25;
figure();
semilogy(MAX_Errors,'.-','MarkerSize',14,'Color','k','linewidth',2);
hold on;
grid on;
box on; 
axis([1 13 10^-17 10^-13]); 
semilogy(MEAN_Errors,'-x','MarkerSize',10,'Color','k','linewidth',2);

set(gca,'FontSize',12) 
set(gca,'Xtick',1:13,'XTickLabel',{'1','','5','','9','','13','','17','','21','','25'});
xlabel('Band-limit, L');
ylabel('E_{max} or E_{mean}');
legend('E_{max}','E_{mean}','Location','northwest');





