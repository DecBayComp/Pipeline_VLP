function [t, diffusivity,  diffusivity_around, diffusivity_variation] = get_evolution_diffusivity_VLP(x_start, y_start, indice_min, indice_max, radius, factor)

t           = zeros(indice_max - indice_min+1,1);
diffusivity = zeros(indice_max - indice_min+1,1);
diffusivity_variation = zeros(indice_max - indice_min+1,1);
diffusivity_around = zeros(indice_max - indice_min+1,1);

for i = indice_min : indice_max

    indice_num = num2str(i);
    load(['Maps_' indice_num '.mat']);
    t(i,1) = Maps(1).t_mean;
    try
    x_centre          = [];
    y_centre          = [];
    for k = 1 : length(Maps);
        x_centre      = [x_centre; Maps(k).center_x];
        y_centre      = [y_centre; Maps(k).center_y];
    end
    
    [indice_voisins, indice_in] = set_neighbors_in_out_boundaries_synapses(x_start, y_start,x_centre, y_centre, factor, radius, Maps);
    [~, indice_in_down]         = set_neighbors_in_out_boundaries_synapses(x_start, y_start,x_centre, y_centre, 1., radius, Maps);
    [indice_voisins_loin,~]     = set_neighbors_in_out_boundaries_synapses(x_start, y_start,x_centre, y_centre, 3., radius, Maps);
    
    indice = 1 ;
    for j =1 : length(indice_in_down)
        diffusivity_in (indice,1) = Maps(indice_in_down(j)).D;
        indice = indice + 1;
    end
    diffusivity_in             = mean(diffusivity_in);
    diffusivity(i,1)            = diffusivity_in;
 
    indice = 1 ;
    for j =1 : length(indice_voisins)
        diffusivity_out(indice,1) = Maps(indice_voisins(j)).D;
        indice = indice + 1;
    end
    diffusivity_out               = mean( diffusivity_out  );
    
    diffusivity_around(i,1)    = diffusivity_out;
    diffusivity_variation(i,1) = (diffusivity_out - diffusivity_in  )./ diffusivity_out; 
    
    clear diffusivity_out diffusivity_loc indice x_centre y_centre indice_voisins;
    clear indice_in indice_in_down  diffusivity_out ;
    clear II r2 x y Maps;
    catch
        diffusivity(i,1)           = nan;
        diffusivity_around(i,1)    = nan;
        diffusivity_variation(i,1) = nan;
    end
end


end