function [t, entropy,  entropy_variation, entropy_all] = get_evolution_entropy_VLP(x_start, y_start, indice_min, indice_max, radius, factor, lambda, sigma, dt_theo)

t                       = zeros(indice_max - indice_min+1,1);
entropy                 = zeros(indice_max - indice_min+1,1);
entropy_variation       = zeros(indice_max - indice_min+1,1);
entropy_all             = zeros(indice_max - indice_min+1,1);

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

            [V, surface ]    = set_extract_potential_global_surface_synapses(Maps,indice_in);
            V                = V - min(V);
            p                = exp(-V)./sum(exp(-V));
            S                = -sum(p.*log(p));
            entropy(i)       = S ;
    %         free_energy(i,1) = F - F_Null;

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
            p         = exp(-V)./sum(exp(-V));
            S_all     = -sum(p.*log(p));

            entropy_all(i,1)       = S_all;
            entropy_variation(i,1) = S_all - S;


            clear indice_voisins indice_in indice_in_down indice_voisins_loin;
            clear x_tot y_tot V K surface S S_all;
            clear x_centre y_centre ; 
        catch
            entropy(i,1)           = nan;
            entropy_all(i,1)       = nan;
            entropy_variation(i,1) = nan;
        end
    else
        entropy(i,1)           = nan;
        entropy_all(i,1)       = nan;
        entropy_variation(i,1) = nan;
    end
end



end