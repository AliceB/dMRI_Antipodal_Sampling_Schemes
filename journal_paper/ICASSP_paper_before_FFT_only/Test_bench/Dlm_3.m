function D_lm = Dlm_3( l,m,m_1,ALPHA, BETA, GAMMA)
% Calculates Wigner D-matrix for a particular zyz rotation defined by (ALPHA, BETA, GAMMA) 

if ( m > l) 
  error('Require m < l.');
end

if ( m_1 > l) 
  error('Require m_1 < l.');
end




D_lm = d_lm( l, m, m_1, BETA).*exp(-1i*m*(ALPHA)).*exp(-1i*m_1*(GAMMA));



end

