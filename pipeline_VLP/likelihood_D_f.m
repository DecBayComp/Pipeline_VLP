function res = likelihood_D_f(dx,dy, dt,param, sigma)
% minus likelihood


D   = param(1);
fx  = param(2);
fy  = param(3);
n   = length(dx);
%sigma = param.sigma;


res = -n*log(4*pi*(D+sigma.^2./dt).*dt) - sum( ((dx-D*fx.*dt).^2 + (dy-D*fy.*dt).^2 ) ./(4*(D+sigma.^2./dt).*dt) ) ;
res = -res;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%