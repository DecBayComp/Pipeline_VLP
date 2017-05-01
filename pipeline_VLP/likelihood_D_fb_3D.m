function res = likelihood_D_fb_3D(dx,dy,dz, dt,param, sigma)
% minus likelihood


D   = param(1);
fx  = param(2);
fy  = param(3);
fz  = param(4);

n   = length(dx);


res = -n*3./2*log(4*pi*(D+sigma.^2./dt)*dt) - sum( ((dx-fx*dt).^2 + (dy-fy*dt).^2 + (dz-fz*dt).^2) ./(4*(D+sigma.^2./dt)*dt) ) ;
res = -res;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%