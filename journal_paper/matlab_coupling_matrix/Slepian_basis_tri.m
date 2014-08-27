function [Basis_mat no_of_basis Basis_eig_value] = Slepian_basis_tri(L,theta_c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Basis_mat = zeros(L^2,L^2);
Basis_eig_value = zeros(1,L^2);
no_of_basis = 0;


for m=0:1:L-1
    

    Gll_mat = Gll_matrix_order_tri(L, m,theta_c );

    [V,D] = eig(Gll_mat);
    D = flipud(diag(D));
    for ii=1:1:length(D)
        %if(abs(D(ii))>0)
           temp_spect = V(:,ii); % choose the temp_spect for which eigenvalue >0.7 
           if m==0
               no_of_basis = no_of_basis + 1;
               Basis_eig_value(no_of_basis) = D(ii);
               for ell=m:1:L-1
                    Basis_mat(ell^2+ell+m+1,no_of_basis) = temp_spect(ell-m+1); 
               end
           else
               no_of_basis = no_of_basis + 1;
               Basis_eig_value(no_of_basis) = D(ii);
               for ell=m:1:L-1
                    Basis_mat(ell^2+ell+m+1,no_of_basis) = temp_spect(ell-m+1); 
               end
               no_of_basis = no_of_basis + 1;
               Basis_eig_value(no_of_basis) = D(ii);
               for ell=m:1:L-1
                    Basis_mat(ell^2+ell-m+1,no_of_basis) = (-1)^m*temp_spect(ell-m+1); 
               end
           end
           
        %end
    end

end

Basis_mat = Basis_mat(:,1:no_of_basis);
Basis_eig_value = Basis_eig_value(1:no_of_basis);

end

