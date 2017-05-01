function [V,surface ] = set_extract_potential_global_surface_synapses(Maps,indice_in)

V = [];
% surface = [];
x_tot = [];
y_tot = [];


for i =1 : length(indice_in)

V     = [V; Maps(indice_in(i)).V];
x_tot = [x_tot; Maps(indice_in(i)).x];
y_tot = [y_tot; Maps(indice_in(i)).y];

end

K         = convhull(x_tot    , y_tot );
surface   = polyarea(x_tot(K) , y_tot(K) );


end