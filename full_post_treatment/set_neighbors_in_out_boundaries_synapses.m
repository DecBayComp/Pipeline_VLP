function [indice_voisins, indice_in] = set_neighbors_in_out_boundaries_synapses(x_start, y_start,x_centre, y_centre, factor, radius, Maps)


theta    = [0:0.1:2*pi];
n_theta  = length(theta);
x_target =  x_start + factor*radius*cos(theta);
y_target =  y_start + factor*radius*sin(theta);
index    = [];



d_centre  = (x_centre - x_start).^2 + (y_centre - y_start ).^2;
indice_in = [];
for i     = 1 : length(Maps)
   
    if (d_centre(i) < (factor*radius).^2 )
        indice_in = [indice_in; i];
    end
    
end

indice_voisins = [];

for i = 1 : n_theta
    
   d             = (x_target(i) - x_centre).^2 + (y_target(i) - y_centre ).^2; 
   [~,I]         = sort(d,'ascend'); 
   indice_voisins = [indice_voisins , I(1:2)]; 
   
end


indice_voisins = unique(indice_voisins);


end
% 
% indice_voisins = [];
% for i = 1 : length(indice_in)
%     indice_voisins =[indice_voisins;  Maps(indice_in(i)).plus_y_index'; ...
%         Maps(indice_in(i)).plus_x_index'; Maps(indice_in(i)).minus_y_index'; ...
%         Maps(indice_in(i)).minus_x_index' ];
% end
% 
% 
% for i = 1 : length(indice_in)
% 
%     II = indice_voisins == indice_in(i);
%     if sum(II) ==0 
%     else
%         indice_voisins = indice_voisins(~II);
%     end
% end
% 
% 
% indice_voisins = unique( indice_voisins);
% for i = 1 : length(indice_voisins)
%     indice_voisins =[indice_voisins;  Maps(indice_voisins(i)).plus_y_index'; ...
%         Maps(indice_voisins(i)).plus_x_index'; Maps(indice_voisins(i)).minus_y_index'; ...
%         Maps(indice_voisins(i)).minus_x_index' ];
% end
% 
% for i = 1 : length(indice_in)
% 
%     II = indice_voisins == indice_in(i);
%     if sum(II) ==0 
%     else
%         indice_voisins = indice_voisins(~II);
%     end
% end
% 
% 
% indice_voisins = unique( indice_voisins);


