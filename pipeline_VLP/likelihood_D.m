function res = likelihood_D(dr2,dt, D, sigma)
% minus likelihood

    n = length(dr2);

    res = -n*log(4*pi*(D+sigma^2./dt)*dt) - sum( dr2./(4*(D+sigma^2./dt)*dt) );
    res = -res;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

