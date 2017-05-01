function [x_target, y_target, intensities] = get_target_densities_VLP(xyt, density_criterion)



sigma = 0.1;
dx    = 0.05;

x = xyt(:,1);
y = xyt(:,2);



[densities, XX, YY] = effective_densities_control(x, y, sigma, dx);
[cent]              = FastPeakFind(densities,density_criterion );

intensities = [];
for j = 1 : floor(length(cent)/2)
    intensities = [intensities; densities(cent(2*j),cent(2*j-1))];
end




indice = 1;
for j = 1 : floor(length(cent)/2)
    if (intensities(j)>= density_criterion)
    x_target(indice,1) =   XX(cent(2*j), cent(2*j-1));
    y_target(indice,1) =   YY(cent(2*j), cent(2*j-1));
    indice             = indice + 1;
    end
end

%plot_local_density_max_VLP( XX, YY,densities,cent);

if ~exist('x_target')
    x_target = [];
    y_target = [];
    
end


end