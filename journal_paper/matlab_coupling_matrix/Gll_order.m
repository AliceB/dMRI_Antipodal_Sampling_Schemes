   function [ Gll ] = Gll_order(ell,ell_1, m,theta_c )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


n_range = abs(ell-ell_1):1:ell+ell_1;
[w3j1]=wigner3jm(ell+ell_1,ell_1,ell,0,0,0);

[w3j2]=wigner3jm(ell+ell_1,ell_1,ell,0,m,-m);
% computes legendre % lazy method; can be made faster

leg_vec = zeros(1,ell+ell_1+3); % for ell=-1, 0, 1...ell+ell_1+1 % 3 extras for -1, 0 and ell+ell+1
leg_vec_diff = zeros(1, ell+ell_1+1);
leg_vec(1) = 1; % for ell=-1

%%
% 1st method
% for ll=0:1:ell+ell_1+1
%    
%     temp = legendre(ll,cos(theta_c));
%     leg_vec(ll+2) = temp(1);
%     
% end


% 2nd method
leg_vec(2) = 1; % for ell=-1
for ll=1:1:ell+ell_1+1
    leg_vec(ll+2) = ((2*ll-1)*cos(theta_c)*leg_vec(ll+1) - (ll-1)*leg_vec(ll))/ll;
    
end
    


%%


for ll=0:1:ell+ell_1
   
    leg_vec_diff(ll+1) = leg_vec(ll+1)- leg_vec(ll+3);
end

prod_temp = w3j1.*w3j2.*leg_vec_diff;
Gll =((-1)^m)*sqrt((2*ell+1)*(2*ell_1+1))/2* sum(prod_temp(abs(ell-ell_1)+1:ell+ell_1+1));
    

end

