function [t, y] = Gaussian_average_heterogeneous_temporal_data_1D(t, x, tau, n_resample)


n     = n_resample;
t_min = min(t);
t_max = max(t);

dt    = (t_max - t_min)./n;
tt    = [t_min: dt:t_max];
tt    = tt';
n_tau = round(tau./dt);

yy    = spline(t,x,tt);
y     = gaussian_smooth(yy, n_tau);


y     = spline(tt,y,t);
end