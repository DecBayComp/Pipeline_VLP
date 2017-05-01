function [V_max , V_min] = set_min_max_V_synapses(Maps ,x_start,y_start,indice_voisins, factor, radius);


for i = 1 : length(Maps)
    
   xx(i,1) = Maps(i).center_x;
   yy(i,1) = Maps(i).center_y;
   VV(i,1) = Maps(i).V; 
    
end

d_loc  = (xx-x_start).^2 + (yy-y_start).^2;
II     = d_loc < (factor.*radius).^2;
VV_loc = VV(II);
V_min  = min(VV_loc);

V_out = [];
for i = 1 : length(indice_voisins)
   V_out = [V_out; Maps(indice_voisins(i)).V];
end
V_max = mean(V_out);




end