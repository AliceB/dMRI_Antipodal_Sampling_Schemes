function [TT_updated min_cond_vec] = nsht_ordered_theta(L)
%need to alter so reads in the .mat file corresponding to each sampling
%scheme

% nsht_ordered_theta - Determine the optimal placement of rings along theta
%
% Either reads in existing .mat files with optimal locations or performs
% step 4 of optimal placement method to create a .mat file with optimal
% locations in it
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
% Original Author: Zubair Khalid
% Altered by: Alice Bates, July 2014
% NSHT package to perform spherical harmonic transforms
%
% for i = 7:4:63
%     nsht_ordered_theta(i);
% end

% Check arguments.
if ~isreal(L)
      error('Harmonic band-limit must be real');
end 


try

    %load(['theta_locations/40L_samples_theta_antipodal_L_',num2str(L),'.mat']);
    load(['theta_locations/V2_just_enough_samples_L_',num2str(L),'.mat']);
    %load('nothing');
    %load(['theta_locations/equator_offset_theta_antipodal_L_',num2str(L),'.mat']);
    return
catch err
    disp('computing optimal ring placement');
%     %optimally place rings 
    
    NUM_RINGS = L; %number of rings to take physical samples over
    TT = nsht_indexed_theta(NUM_RINGS);
    TT_ordered_index = [NUM_RINGS];
   
    %check that legendre polynomial of degree L-1 and order L-2 not zero
    %for this sample
    [P_mat Sc] = nsht_legmat_mex(TT(TT_ordered_index), L, L-2);
    PP = 10.^Sc.*P_mat;
    PP = PP(2);
    
    i=1;
    while (abs(PP) < 1e-14)
        disp('PP SMALL --------------------------------------');
        TT_ordered_index = [NUM_RINGS - i];
        [P_mat Sc] = nsht_legmat_mex(TT(TT_ordered_index), L, L-2);
        PP = 10.^Sc.*P_mat;
        PP = PP(2);
        i = i+1;
    end
    
    min_cond_vec = [];

    for m = L-3:-2:0
            index_1 =1:2:(L-m);
            index_2 = 2:2:L-(m-1);
            
            [P_mat_1 Sc_1] = nsht_legmat_mex(TT, L, m); 
            
            Y_mat_1 = zeros(length(index_1),length(index_1));
            

            PP_1 = 10.^Sc_1.*P_mat_1;
            PP_1 = PP_1(index_1,:); %keep only even degree rows

            
            for jj=1:1:length(TT_ordered_index)
                Y_mat_1(:,jj) = PP_1(:,TT_ordered_index(jj));
            end

            
             
            if m~=0  % work out P matrix for m-1
                [P_mat_2 Sc_2] = nsht_legmat_mex(TT, L, m-1); %m odd            
                Y_mat_2 = zeros(length(index_2),length(index_2));
             
                PP_2 = 10.^Sc_2.*P_mat_2;
                PP_2 = PP_2(index_2,:); %keep only even degree rows
               for jj=1:1:length(TT_ordered_index)
                    Y_mat_2(:,jj) = PP_2(:,TT_ordered_index(jj));
                end
             end
            

            ii_vec=[];
            cond_vec_1 = [];
            cond_vec_2 = [];
            for ii=1:1:length(TT)

                if isempty(find(TT_ordered_index==ii))
                    Y_mat_1(:,end) = PP_1(:,ii);
                    ii_vec = [ii_vec ii];
                    cond_vec_1 = [cond_vec_1 cond(Y_mat_1)];
                     if m~=0
                         Y_mat_2(:,end) = PP_2(:,ii);
                        cond_vec_2 = [cond_vec_2 cond(Y_mat_2)];
                     end
                end
            end
            
            
            if m~=0
              % Minimising sum of two condition numbers
               [M,II] =min(cond_vec_1+ cond_vec_2);

               
            else
                [M,II] =min(cond_vec_1);
            end
            
            min_cond_vec = [min_cond_vec M];
            TT_ordered_index = [ ii_vec(II) TT_ordered_index];

    end

    TT_updated = TT(TT_ordered_index);
    
    save(['theta_locations\V2_just_enough_samples_L_' num2str(L) '.mat'], 'TT_updated', 'min_cond_vec'); 
   
end %% end try/catch

end


