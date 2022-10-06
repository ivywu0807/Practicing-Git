function [f ,F] = asymstab(xvec,a,b)
% pdf and , optionally , cdf of the asymmetric stable . See also asymstabpdf .m % Nolan's routine: Call stablepdf(xvec,[a,b],1)
if nargin<3, b=0; end
bordertol=1e-8; lo=bordertol ; hi=1-bordertol ; tol=1e-7; xl=length(xvec); F=zeros(xl,1); f=F;
for loop=1:length(xvec)
    x=xvec(loop); f(loop)= quadl(@fff,lo,hi,tol,[],x,a,b,1) / pi; 
    if nargout>1
        F(loop)=0.5-(1/pi)* quadl(@fff,lo,hi,tol,[],x,a,b,0); end
end;
function I = fff (uvec,x,a,b,dopdf) 
subs = 1; I = zeros(size(uvec)); 
for ii=1:length(uvec)
    u=uvec ( ii ) ;
    if subs==1, t =(1-u)/u; else t =u/(1-u); end
    if a==1, cf = exp( -abs(t)*( 1 + 1i*b*(2/pi)*sign(t)*log(t)) ); 
    else cf = exp( -((abs(t))^a) *( 1 - 1i*b*sign(t)*tan(pi*a/2)) );
    end
    z = exp(-1i*t*x).* cf; if dopdf==1, g=real(z); else g=imag(z)./t; end
    if subs==1, I( ii )=g*u^(-2); else I( ii )=g*(1-u)^(-2); end 
end

