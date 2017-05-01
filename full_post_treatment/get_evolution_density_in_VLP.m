function [t, density] = get_evolution_density_in_VLP(x_start, y_start, indice_min, indice_max, radius, factor)

t       = zeros(indice_max - indice_min+1,1);
density = zeros(indice_max - indice_min+1,1);

for i = indice_min : indice_max

    indice_num = num2str(i);
    load(['Maps_' indice_num '.mat']);
    t(i,1) = Maps(1).t_mean;
    
    x = [];
    y = [];
    for  j = 1 : length(Maps)
       x = [x; Maps(j).x];
       y = [y; Maps(j).y];
    end
    
    r2 = (x-x_start).*(x-x_start) + (y-y_start).*(y-y_start);
    II = r2 <= factor*factor*radius*radius;
    
    density(i,1) = sum(II)./(pi*radius*radius);
    clear II r2 x y Maps;
end


end