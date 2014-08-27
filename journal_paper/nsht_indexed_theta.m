function [ TT ] = nsht_indexed_theta(L)
%#codegen
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

L=L-1; % convention - bandlimit
TT_temp = pi*(2*(0:1:L)+1)/(2*L+1);



% %%%%%%%%%%%%%%%%%%%%%%%%
% TT_temp = pi*(2*(0:1:L)+1)/(2*L+1);
% TT_new = zeros(1,111);
% TT_new(end) = pi;
% 
% for i=0:1:length(TT_new)-2
% TT_new(end-i-1) = TT_new(end-i) - sin(TT_temp(i+1));
% end
% 
% TT_temp = (pi*(TT_new-min(TT_new))/(max(TT_new)-min(TT_new)));
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


TT = zeros(size(TT_temp));

index_middle = ceil((L+1)/2);

TT(1) = TT_temp(index_middle);
%% L even, number of samples are odd. After first one should go left
if (mod(L,2)==0)
    for i=index_middle+1:1:L+1
     TT(2*(i-index_middle)+1) = TT_temp(i);
    end
    
    for i=index_middle-1:-1:1

     TT(2*(index_middle-i)) = TT_temp(i);
    end
    
else
% L odd, number of samples are even. After first one should go right    
    for i=index_middle+1:1:L+1

     TT(2*(i-index_middle)) = TT_temp(i);
    end
    
    for i=index_middle-1:-1:1

     TT(2*(index_middle-i)+1) = TT_temp(i);
    end

end

TT = fliplr(TT);

end

