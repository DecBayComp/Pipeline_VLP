function [x_start, y_start, II, JJ] = give_center_VLP(Maps, value)



[ ~, ~, densities, XX, YY, ~, ~] = set_seed_densities_VLP(Maps);
[cent]                           = FastPeakFind(densities,value );

if isempty(cent)
   x_start = [];
   y_start = [];
   II      = [];
   JJ      = [];
else
   II      = cent(1:2:end);
   JJ      = cent(2:2:end);
   n_II    = length(II);
   x_start = zeros(n_II,1);
   y_start = zeros(n_II,1);
   for i = 1 : n_II
    x_start(i,1) = XX(JJ(i),II(i));
    y_start(i,1) = YY(JJ(i),II(i));
   end
    
end










end