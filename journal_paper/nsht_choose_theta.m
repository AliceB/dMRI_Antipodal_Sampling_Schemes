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

for L = 9:4:61
    if mod(L,2) ==0
        L = L -1;
    end
   
        TT_temp = pi*(0.00000001:(THETA_FAC*L))/((THETA_FAC*L)) ; 
        LL = length(TT_temp);
        
        %choose first sample at north pole
        index_final = 1;
        
        for ell = 3:2:L
            m=0;
            [P, Sc] = nsht_legmat_mex(TT_temp, ell, m);
            PP = 10.^Sc.*P;
            PP = PP(1:2:end,:); %get rid of odd legendre polyomials 
            cond_vec = 10^18*ones(1,LL);

            for ii=1:1:LL
                  cond_temp = cond( transpose(PP(:,[ii, index_final])));
                    if cond_temp < min(cond_vec)
                        c_temp = cond_temp;
                        II_final = ii;

                    end
                    cond_vec(ii) = cond_temp; 

            end

            index_final = [II_final, index_final]; 

        end

        %%
        P_mat = transpose(PP(:,index_final));
        cond(P_mat)

        %%
        theta_final = TT_temp(index_final);
        %%save the theta locations
        str = ['theta_locations/P0condmin_theta_dMRI_Journal_',num2str(L),'.mat'];
        save(str,'theta_final'); 
  
end