function [Maps_densities] = output_density_from_Maps(Maps, sigma,dx)

if isfield(Maps, 'x_all')
    
    x = Maps(1).x_all;
    y = Maps(1).y_all;
    
else

    x = [];
    y = [];
    for i = 1 : length(Maps)
        x = [x;Maps(i).x];
        y = [y;Maps(i).y];
    end
end



[densities, XX, YY] = generate_density_from_point_map(x, y, sigma,dx);

Maps_densities.XX        = XX;
Maps_densities.YY        = YY;
Maps_densities.densities = densities;







