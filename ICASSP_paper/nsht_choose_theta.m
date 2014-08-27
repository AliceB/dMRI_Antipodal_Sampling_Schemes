%steps1_to_3_sampling_scheme_design.m carries out steps 1 to 3 of the optimal 
%placement method contained in the paper Novel Sampling Scheme on the Sphere for
%Head-Related Transfer Function Measurements.   
%This trial starts with one sample position and computes the next sample position which
%minimizes the condition number of 2x2 matrix, choosen next which
%minimizes 3x3 and so on.
%A .mat file containing theta locations is created 
%author: Alice Bates
%date: July 2014

%% parameters
THETA_FAC = 90; %how many spare rings to choose from in multiplicative terms

for L = 9 %8:4:64
   
        TT_temp = pi*(0.00000001:(THETA_FAC*L))/((THETA_FAC*L)) ; 
        TT_temp = [TT_temp, (pi - TT_temp)];

        LL = length(TT_temp)/2;
        
        % L odd
        if mod(L,2) == 1 
            %choose first sample at north pole
            index_final = 1;
            ell_values = 3:2:L;
            %%
        else %L even 
             ell_values = 2:2:L;
            index_final = [];
        end
 
        for ell = ell_values
            m=0;
            [P, Sc] = nsht_legmat_mex(TT_temp, ell, m);
            PP = 10.^Sc.*P;

            cond_vec = 10^18*ones(1,LL);

            for ii=1:1:LL
                  cond_temp = cond( transpose(PP(:,[ii + LL,ii, index_final])));
                    if cond_temp < min(cond_vec)
                        c_temp = cond_temp;
                        II_final = ii;

                    end
                    cond_vec(ii) = cond_temp; 

            end

            index_final = [(II_final +LL), II_final, index_final]; 

        end

        %%
        P_mat = transpose(PP(:,index_final));
        cond(P_mat)

        %%
        theta_final = TT_temp(index_final);
        %%save the theta locations
        str = ['../Data/P0condmin_theta_dMRI_ICASSP_',num2str(L),'.mat'];
        save(str,'theta_final'); 
  
end