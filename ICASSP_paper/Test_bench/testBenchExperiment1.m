%testBench.m generates reconstruction error plots for antipodal signal
%author: Alice Bates
%date: July 2014

% Install SHT 
install_sht

MAX_Errors = zeros(1,15);
MEAN_Errors = zeros(1,15);

% 
global conditionNum ;
conditionNum = zeros(1,L);


index = 1; 
for L = 8:4:64
    for j =1:10 %repeat experiment 10 times for each L

        %need to consider case L odd
        %% Generate random antipodal signal from L*(L+1)/2 samples
%         ft = zeros(1,L^2);
%         for l = L-1:-2:0
%              ft(l^2+1:(l+1)^2) =  rand(1,(2*l+1)) + 1i*rand(1,(2*l+1));
%         end    
%         
%         ft_temp = ft;
        
        %copy samples to antipodal locations
        [~, phis] = nsht_sampling_points(L);
        for l = L-2:-2:0
            %search through ring antipodal to ring l (next largest ring)
           for j = l^2+1:(l+1)^2
                for i = l^2+2*l+2:(l+2)^2
                   %if antipodal point
                  if phis(j) == mod((phis(i)+pi),2*pi)
                     ft(j) = ft(i); 
                  end
                end
           end
        end
       
        %% forward transform
        flmt = nsht_forward(ft,L);

        %% inverse transform
        fr  = nsht_inverse(flmt,L);


        %% Absolute error
        Error_max = max(abs(ft-fr));

        %%MEAN error
        Error_mean = sum(abs(ft-fr))/L^2;  
        
        %store average RMS and MAX error for plotting
        MAX_Errors(index) = MAX_Errors(index) +  Error_max*0.1;
        MEAN_Errors(index) = MEAN_Errors(index) + Error_mean*0.1;

   end
    
    index = index + 1;
end



%% plot max and mean error
x = 8:4:64;
figure();
semilogy(MAX_Errors,'.-','MarkerSize',14,'Color','k','linewidth',2);
hold on;
semilogy(MEAN_Errors,'-x','MarkerSize',10,'Color','k','linewidth',2);

set(gca,'FontSize',12) 
set(gca,'Xtick',1:15,'XTickLabel',{'8','','16','','24','','32','','40','','48','','56','','64'});
xlabel('Band-limit, L');
ylabel('E_{max} or E_{mean}');
legend('E_{max}','E_{mean}','Location','northwest');
