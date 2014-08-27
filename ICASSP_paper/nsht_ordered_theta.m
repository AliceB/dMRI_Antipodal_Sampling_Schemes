function [ TT_updated min_cond_vec] = nsht_ordered_theta(L)
% nsht_ordered_theta - Determine the optimal placement of rings along theta
%
% 
%
% Default usage is given by
%
%   f = nsht_ordered_theta(L)
%
% where L is the harmonic band-limit. TT_updated denotes the \theta vector 
% contains sampling positions along \theta. min_cond_vector contains the condition number
% of the matrix Y_m. See the following paper for details.
%
%
% Author: Zubair Khalid
%
% NSHT package to compute spherical harmonic transform of band-limited
% signal
% Copyright (C) 2014  Zubair Khalid
% See LICENSE.txt for license details
%%
% Check arguments.
if ~isreal(L)
      error('Harmonic band-limit must be real');
end


try
    
   load(['../Data/theta_minPm_dMRI_ICASSP_' num2str(L)], 'TT_updated', 'min_cond_vec');
  %load(['../Data/bothOptim_theta_minPm_dMRI_ICASSP_' num2str(L)], 'TT_updated', 'min_cond_vec');
 
  
    return
catch err
%     %THETA_FAC = 40; %determines the number of rings have to choose from 
%     disp('computing optimal ring placement');
%     %optimally place rings 
%     TT = nsht_indexed_theta(THETA_FAC*L);
%     TT_ordered_index = [THETA_FAC*L, 2*THETA_FAC*L]; %choose the ring closest to the equator
%   
%    min_cond_vec = [];
%     %%
% 
%     for m = L-3:-2:0
%             
%             [P1, Sc1] = nsht_legmat_mex(TT, L, m);     
%             Y_mat_1 = zeros(L-m,L-m);
%             PP_1 = 10.^Sc1.*P1;
% 
%             for jj=1:1:length(TT_ordered_index)
%                 Y_mat_1(:,jj) = PP_1(:,TT_ordered_index(jj));               
%             end
% 
%             if m~=0
%                 [P2, Sc2] = nsht_legmat_mex(TT, L, m-1); 
%                 Y_mat_2 = zeros(L-(m-1),L-(m-1));
%                 PP_2 = 10.^Sc2.*P2;
%                 for jj=1:1:length(TT_ordered_index)
%                     Y_mat_2(:,jj) = PP_2(:,TT_ordered_index(jj));
%                 end
%                 cond_vec_2 = [];
%             end
%             ii_vec=[];
%             cond_vec_1 = [];
%             
%             
%             for ii=1:1:length(TT)/2
%                 if isempty(find(TT_ordered_index==ii))
%                     Y_mat_1(:,end) = PP_1(:,ii);
%                     ii_vec = [ii_vec ii];
%                     cond_vec_1 = [cond_vec_1 cond(Y_mat_1)];
%                     if m~=0
%                         Y_mat_2(:,end-1)= PP_2(:,ii);
%                         Y_mat_2(:,end) = PP_2(:,ii + THETA_FAC*L);
%                         cond_vec_2 = [cond_vec_2 cond(Y_mat_2)];
%                      end
%                         
%                 end
%             end
%             
%             if m~=0
%                 [M,II] = min(cond_vec_1+ cond_vec_2);
%             else
%                 [M,II] =min(cond_vec_1);
%             end
%             min_cond_vec = [min_cond_vec M];
%             TT_ordered_index = [ (ii_vec(II) + THETA_FAC*L) ii_vec(II) TT_ordered_index];      
% 
%     end
% 
%     TT_updated = TT(TT_ordered_index);
%     %%save theta locations
%     save(['../Data/theta_minPm_dMRI_ICASSP_' num2str(L) '.mat'], 'TT_updated', 'min_cond_vec');
% end %% end try/catch


%% change location of rings such as to minimise condition number of Pm larger m
    load(['../Data/P0condmin_theta_dMRI_ICASSP_',num2str(L),'.mat']);
    TT = theta_final;  
   [~,temp_index ] = min(abs(theta_final - 0.5*pi)); %find location of the largest ring (closest to pi/2)
   if mod(temp_index,2) == 0
        TT_ordered_index = [temp_index-1 temp_index];
   else
        TT_ordered_index = [temp_index+1 temp_index];
   end
   
   min_cond_vec = [];
    %%

    for m = L-3:-2:0
            
        [P1, Sc1] = nsht_legmat_mex(TT, L, m);     
        Y_mat_1 = zeros(L-m,L-m);
        PP_1 = 10.^Sc1.*P1;

        for jj=1:1:length(TT_ordered_index)
            Y_mat_1(:,jj) = PP_1(:,TT_ordered_index(jj));               
        end

        if m~=0
            [P2, Sc2] = nsht_legmat_mex(TT, L, m-1); 
            Y_mat_2 = zeros(L-(m-1),L-(m-1));
            PP_2 = 10.^Sc2.*P2;
            for jj=1:1:length(TT_ordered_index)
                Y_mat_2(:,jj) = PP_2(:,TT_ordered_index(jj));
            end
            cond_vec_2 = [];
        end
        ii_vec=[];
        cond_vec_1 = [];


        for ii=1:1:length(TT)
            if isempty(find(TT_ordered_index==ii))
                Y_mat_1(:,end) = PP_1(:,ii);
                ii_vec = [ii_vec ii];
                cond_vec_1 = [cond_vec_1 cond(Y_mat_1)];
                if m~=0
                    Y_mat_2(:,end-1)= PP_2(:,ii);
                    if mod(ii,2) == 1
                        Y_mat_2(:,end) = PP_2(:,ii+1);
                    else
                        Y_mat_2(:,end) = PP_2(:,ii-1);
                    end
                    cond_vec_2 = [cond_vec_2 cond(Y_mat_2)];
                 end

            end
        end

        if m~=0
            [M,II] = min(cond_vec_1+ cond_vec_2);
        else
            [M,II] =min(cond_vec_1);
        end
        min_cond_vec = [min_cond_vec M];
        if mod(ii_vec(II),2) == 1
            TT_ordered_index = [ (ii_vec(II) + 1), ii_vec(II), TT_ordered_index];      
        else
            TT_ordered_index = [ (ii_vec(II) - 1), ii_vec(II), TT_ordered_index]; 
        end
    end

    TT_updated = TT(TT_ordered_index);
    %%save theta locations
    save(['../Data/bothOptim_theta_minPm_dMRI_ICASSP_' num2str(L) '.mat'], 'TT_updated', 'min_cond_vec');
end %% end try/catch

