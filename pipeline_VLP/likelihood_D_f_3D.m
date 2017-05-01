function res = likelihood_D_f_3D(dx,dy,dz, dt,param, sigma)
% minus likelihood


D   = param(1);
fx  = param(2);
fy  = param(3);
fz  = param(4);

n   = length(dx);
%sigma = param.sigma;


res = -n.*3./2*log(4*pi*(D+sigma.^2./dt).*dt) - sum( ((dx-D*fx.*dt).^2 + (dy-D*fy.*dt).^2 + (dz-D*fz.*dt).^2) ./(4*(D+sigma.^2./dt).*dt) ) ;
res = -res;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%