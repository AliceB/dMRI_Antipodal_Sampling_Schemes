function fact_c = norm_factor(l,m)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

L1 = l+m;
L2 = l-m+1;

num = sqrt((2*l+1)/(4*pi));
den = 1/num;

for i= L2:1:L1
    
    den = den *sqrt(i);
    

end


fact_c = 1./den;
end

