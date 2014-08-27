function [ Gll_mat ] = Gll_matrix_order_tri(L, m,theta_c )

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



Gll_mat = zeros(L-m,L-m);

% tri-diagonal terms only EXCEPT FOR the last enrty
for ell=m:1:L-2
   
  
    Gll_mat(ell+1-m,ell+1-m)=-ell*(ell+1)*cos(theta_c);
    temp = (ell*(ell+2) - (L-1)*(L+1))*sqrt(((ell+1)^2-m^2)/((2*ell+1)*(2*ell+3)));
    Gll_mat(ell+2-m,ell+1-m)=temp;
    Gll_mat(ell+1-m,ell+2-m)=temp;

end
   
ell = L-1;
Gll_mat(ell+1-m,ell+1-m)=-ell*(ell+1)*cos(theta_c);



end