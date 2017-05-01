function [t, potential_depth,  potential_depth_far, potential_depth_global] = get_evolution_potential_depth_VLP(x_start, y_start, indice_min, indice_max, radius, factor, lambda, sigma, dt_theo)

t                       = zeros(indice_max - indice_min+1,1);
potential_depth         = zeros(indice_max - indice_min+1,1);
potential_depth_far     = zeros(indice_max - indice_min+1,1);
potential_depth_global  = zeros(indice_max - indice_min+1,1);


for i = indice_min : indice_max
    fprintf('%i\n', i);
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
    %         [~, indice_in_down]         = set_neighbors_in_out_boundaries_synapses(x_start, y_start,x_centre, y_centre, 1., radius, Maps);
            [indice_voisins_loin,~]  = set_neighbors_in_out_boundaries_synapses(x_start, y_start,x_centre, y_centre, 3., radius, Maps);


            [V_max , V_min]          = set_min_max_V_synapses(Maps ,x_start,y_start,indice_voisins, factor, radius);
            dV                       = V_max - V_min;
            potential_depth(i,1)     = dV;

            [V_max , V_min]          = set_min_max_V_synapses(Maps ,x_start,y_start,indice_voisins_loin, factor, radius);
            dV_far                   = V_max - V_min;
            potential_depth_far(i,1) = dV_far;


            V     = [];
            x_tot = [];
            y_tot = [];
            for j = 1: length(Maps)
                V(j)  = Maps(j).V ;
                x_tot = [x_tot; Maps(j).x];
                y_tot = [y_tot; Maps(j).y];
            end

            dV_global                   = mean(V) - V_min;
            potential_depth_global(i,1) = dV_global;


            clear indice_voisins indice_in indice_in_down indice_voisins_loin;
            clear x_tot y_tot V K surface F F_all;
            clear x_centre y_centre ; 
        catch
            potential_depth(i,1)           = nan;
            potential_depth_far(i,1)       = nan;
            potential_depth_global(i,1)    = nan;
            
        end
    else
        potential_depth(i,1)           = nan;
        potential_depth_far(i,1)       = nan;
        potential_depth_global(i,1)    = nan;
    end
end


end