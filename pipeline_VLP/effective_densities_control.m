function [densities, XX, YY] = effective_densities_control(x, y, sigma, dx, varargin);

%% check if multiple processors can be used
% selective_kill_parallel_multi_version;
% start_check_parallel_multi_version;


if nargin <3
    sigma = 0.1;
end
sigma2 = sigma.^2;

% dx = 0.025;
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

for i = 1 : length(x)

    r2 = (XX - x(i)).^2 + (YY - y(i)).^2;
    densities = densities + exp(-r2./(2*sigma2))./(2*pi*sigma2);
    
end

%densities = densities*dx^2;
