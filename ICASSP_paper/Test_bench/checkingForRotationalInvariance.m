%testBench.m generates reconstruction error plots for antipodal signal used
%in ICASSP paper
%author: Alice Bates
%date: September 2014

% Install SHT 
install_sht
NUM_BAND_LIMITS = 13;
MAX_BAND_LIMIT = 25;
NUM_ROTATIONS = 6;
MAX_Errors = zeros(NUM_ROTATIONS,NUM_BAND_LIMITS);
MEAN_Errors = zeros(NUM_ROTATIONS,NUM_BAND_LIMITS);
NUM_EXPERIMENTS = 10;

%generate random rotations
alpha = 2*pi*rand(1,NUM_ROTATIONS-1);
alpha = [alpha,0];
beta = pi*rand(1,NUM_ROTATIONS-1);
beta = [beta, 0];
gamma  = 2*pi*rand(1,NUM_ROTATIONS-1);
gamma = [gamma, 0];
rotation_matrices = zeros(MAX_BAND_LIMIT,MAX_BAND_LIMIT,NUM_ROTATIONS);


%create rotation matrices for largest band-limit
for k=1:NUM_ROTATIONS
    rowIndex = 1;
    colIndex = 1;
    for l=0:MAX_BAND_LIMIT
        for m=-l:l
            for m_1 = -l:l
                rotation_matrices(rowIndex,colIndex+l+m_1,k) = Dlm_3( l,m,m_1,alpha(k), beta(k), gamma(k));
            end
            rowIndex = rowIndex + 1;
        end
        colIndex = colIndex + 2*l +1;
    end
end
    
index = 1; 

for L = 1:2:25
    for j =1:NUM_EXPERIMENTS %repeat experiment x times for each L


        %% Generate random antipodal signal
        flmt = zeros(1,L^2);
        for l = 1:2:L
             flmt(l^2-2*l+2:l^2) =  rand(1,(2*l-1)) + 1i*rand(1,(2*l-1));
        end
        
        for k=1:NUM_ROTATIONS
            
            %work out coefficients of rotated signal
            flmt_rot = rotation_matrices(1:L^2,1:L^2,k)*flmt.';
            flmt_rot = flmt_rot.';
    
            %% inverse transform
            ft_rot  = nsht_inverse(flmt_rot,L);


            %% forward transform
            flmr_rot = nsht_forward(ft_rot,L);


            %% Absolute error
            Error_max = max(abs(flmt_rot-flmr_rot));

            %%MEAN error
            Error_mean = sum(abs(flmt_rot-flmr_rot))/L^2;

            %store average RMS and MAX error for plotting
            MAX_Errors(k,index) = MAX_Errors(k,index) +  Error_max*(1/NUM_EXPERIMENTS);
            MEAN_Errors(k,index) = MEAN_Errors(k,index) + Error_mean*(1/NUM_EXPERIMENTS);
        end

   end
    
    index = index + 1;
end
%%
%legend strings
str1 = strcat('(',num2str(alpha(1),2),',', num2str(beta(1),2),',' ,num2str(gamma(1),2),')');
str2 = strcat('(',num2str(alpha(2),2),',', num2str(beta(2),2),',' ,num2str(gamma(2),2),')');
str3 = strcat('(',num2str(alpha(3),2),',', num2str(beta(3),2),',' ,num2str(gamma(3),2),')');
str4 = strcat('(',num2str(alpha(4),2),',', num2str(beta(4),2),',' ,num2str(gamma(4),2),')');
str5 = strcat('(',num2str(alpha(5),2),',', num2str(beta(5),2),',' ,num2str(gamma(5),2),')');
str6 = strcat('(',num2str(alpha(6),2),',', num2str(beta(6),2),',' ,num2str(gamma(6),2),')');
%% plot mean error
x = 1:2:25;
figure();

 

%colours = [1 0 0;0 1 0;1 1 0;0 0 0 ;0 1 1;1 0 1;0 0 1;0.5 0.25 0;0 0.5 0.75 ; 0.5 0.9 0.1];
%colours =[1 0 0;0 1 0;0 0 1; 0 1 1 ; 1 0 1];
semilogy(MEAN_Errors(1,:),'-x','MarkerSize',10,'Color','k','linewidth',2);

hold on;axis([1 13 10^-17 10^-14]); 
grid on;
box on;
semilogy(MEAN_Errors(2,:),'-o','MarkerSize',10,'Color','k','linewidth',2);
semilogy(MEAN_Errors(3,:),'-*','MarkerSize',10,'Color','k','linewidth',2);
semilogy(MEAN_Errors(4,:),'-v','MarkerSize',10,'Color','k','linewidth',2);
semilogy(MEAN_Errors(5,:),'-s','MarkerSize',10,'Color','k','linewidth',2);
semilogy(MEAN_Errors(6,:),'-.','MarkerSize',10,'Color','r','linewidth',2);
set(gca,'FontSize',12) 
set(gca,'Xtick',1:NUM_BAND_LIMITS,'XTickLabel',{'1','','5','','9','','13','','17','','21','','25'});
xlabel('Band-limit, L');
ylabel('E_{mean}');
legend(str1,str2,str3,str4,str5,str6,'Location','northwest');


% x = 1:2:25;
% figure();
% hold on;
% grid on;
% box on; 
% %axis([1 13 10^-17 10^-14]); 
% colours = [1 0 0;0 1 0;1 1 0;0 0 0 ;0 1 1;1 0 1;0 0 1;0.5 0.25 0;0 0.5 0.75 ; 0.5 0.9 0.1];
% for k=1:NUM_ROTATIONS
%     semilogy(MAX_Errors(k,:),'.-','MarkerSize',14,'Color',colours(k,:),'linewidth',2);
%    
% end
% set(gca,'FontSize',12) 
% set(gca,'Xtick',1:NUM_BAND_LIMITS,'XTickLabel',{'1','','5','','9','','13','','17','','21','','25'});
% xlabel('Band-limit, L');
% ylabel('E_{max}');
% legend(str1,str2,str3,str4,str5,'Location','northwest');
% 


