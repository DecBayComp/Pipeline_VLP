function res = likelihood_D_fb(dx,dy, dt,param, sigma)
% minus likelihood


D   = param(1);
fx  = param(2);
fy  = param(3);
n   = length(dx);
dt  = dt*ones(n,1);


res = -n*log(4*pi*(D+sigma.^2./dt(1)).*dt(1)) - sum( ((dx-fx*dt).^2 + (dy-fy*dt).^2 ) ./(4*(D+sigma.^2./dt).*dt) ) ;
res = -res;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%