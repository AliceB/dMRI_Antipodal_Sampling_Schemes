function [ Gll_mat ] = Gll_matrix_order(L, m,theta_c )

% computes coupling matrix for the fixed order.

Gll_mat = zeros(L-m,L-m);


for ell=m:1:L-1
    
   for ell_1=ell:1:L-1
   
    Gll = Gll_order(ell,ell_1,m,theta_c);
    Gll_mat(ell+1-m,ell_1+1-m)=Gll;
    Gll_mat(ell_1+1-m,ell+1-m)=Gll;
   end
end
    
    


end

