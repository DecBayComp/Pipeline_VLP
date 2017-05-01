function [t, free_energy,  free_energy_all, free_energy_variation] = get_evolution_free_energy_VLP(x_start, y_start, indice_min, indice_max, radius, factor, lambda, sigma, dt_theo)

t                       = zeros(indice_max - indice_min+1,1);
free_energy             = zeros(indice_max - indice_min+1,1);
free_energy_variation   = zeros(indice_max - indice_min+1,1);
free_energy_all         = zeros(indice_max - indice_min+1,1);

for i = indice_min : indice_max

    indice_num = num2str(i);
    load(['Maps_' indice_num '.mat']);
    Maps = Set_potential_stochastic_integral(Maps, lambda, sigma, dt_theo);
    t(i,1) = Maps(1).t_mean;
    if isfield(Maps,'V')
        
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

            [V, surface ]    =  set_extract_potential_global_surface_synapses(Maps,indice_in);
             V               = V - min(V);
             F               = -log( sum(exp(-V)) )- ( log(surface)) ;
             F_Null          = -log( sum(exp(-V*0.)) )- ( log(surface)) ;
            free_energy(i,1) = F - F_Null;

            V     = [];
            x_tot = [];
            y_tot = [];
            for j = 1: length(Maps)
                V(j)  = Maps(j).V ;
                x_tot = [x_tot; Maps(j).x];
                y_tot = [y_tot; Maps(j).y];
            end
            K         = convhull(x_tot    , y_tot );
            surface   = polyarea(x_tot(K) , y_tot(K) );
            V         = V - min(V);
            F_all     = -log( sum(exp(-V)) ) - log(surface);

            free_energy_all(i,1)       = F_all;
            free_energy_variation(i,1) = F_all - F;


            clear indice_voisins indice_in indice_in_down indice_voisins_loin;
            clear x_tot y_tot V K surface F F_all;
            clear x_centre y_centre ; 
        catch
            free_energy(i,1)           = nan;
            free_energy_all(i,1)       = nan;
            free_energy_variation(i,1) = nan; 
        end
    else
        free_energy(i,1)           = nan;
        free_energy_all(i,1)       = nan;
        free_energy_variation(i,1) = nan;
    end
end


end