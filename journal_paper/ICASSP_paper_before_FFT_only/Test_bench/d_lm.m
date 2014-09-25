function [ dlm ] = d_lm( l, m, m_1, b)
%UNTITLED Summary of this function goes here
%   Only b (angle in radians) can be a vector 

%choosing lower limit of j s.t. denominator positive
if (m>=m_1)
    
   j=0; 
    
else
   j=m_1-m; 
    
end

result =0;


den1 = (l+m_1-j);

den2 = (l-m-j);

den3 = (j-m_1+m);


    while(den1>=0 && den2>=0 && den3>=0 )
                
        
            num = sqrt(factorial(l+m_1)*factorial(l-m_1)*factorial(l+m)*factorial(l-m));

            den = factorial(den1)*factorial(den2)*factorial(den3)*factorial(j);
            
            term1 = (-1)^den3;
            
            term2 = num./den;
            
            term3 = (cos(b./2)).^(2*l-2*j+m_1-m);

            term4 = (sin(b./2)).^(2*j-m_1+m);

            result = result + term1.*term2.*term3.*term4;

            j = j+1;




        


        den1 = (l+m_1-j);

        den2 = (l-m-j);

        den3 = (j-m_1+m);

        
    end


        dlm = result;
        
end

