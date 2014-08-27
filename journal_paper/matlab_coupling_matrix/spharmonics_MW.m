
function [Ylm]=spharmonics_MW(l,m,L,show)

% Define constants (REQUIRED THAT L(DEGREE)>=M(ORDER))
if nargin==0
  l=3;   % DEGREE
  m=2;   % ORDER
  L=100;
end

if nargin<3
  L=100;
end

if l<m, error('The ORDER (m) must be less than or eqaul to the DEGREE(l).'); end
%%

TT = pi*(2*(0:1:L)+1)/(2*L+1);

FF = pi*(2*(0:1:2*L))/(2*L+1);


leg_all = legendre(l,cos(TT));


%%

Ylm=norm_factor(l,abs(m)).*transpose(leg_all(abs(m)+1,:))*exp(1i*m*FF);

if m<0
    
   Ylm = (-1)^m*Ylm; 
end




if show
    [FI,THETA]=meshgrid(FF,TT); % creates grid

    Xm = 1*sin(THETA).*cos(FI);
Ym = 1*sin(THETA).*sin(FI);
Zm = 1*cos(THETA);

figure;
h = surf(Xm, Ym, Zm, real(Ylm));
set(h, 'LineStyle','none');


figure;
h = surf(Xm, Ym, Zm, imag(Ylm));
set(h, 'LineStyle','none');

end

 return