function [densities, XX, YY] = effective_densities(x, y, sigma, varargin);

%% check if multiple processors can be used
%     if (matlabpool('size') == 0)
%     myCluster = parcluster('local');
%     matlabpool(myCluster, myCluster.NumWorkers);
%     end


if nargin <3
    sigma = 0.1;
end
sigma2 = sigma.^2;

dx = 0.025;
xmin = min(x);
ymin = min(y);
xmax = max(x);
ymax = max(y);


xx = [xmin:dx:xmax];
yy = [ymin:dx:ymax];


[XX,YY] = meshgrid(xx,yy);


densities = XX*0.;

n_xx = length(xx);
n_yy = length(yy);
%parfor i = 1 : length(x)
for i = 1 : length(x)
     
    
    r2 = (XX - x(i)).^2 + (YY - y(i)).^2;
    densities = densities + exp(-r2./(2*sigma2))./(2*pi*sigma2);
    
end

%densities = densities*dx^2;
