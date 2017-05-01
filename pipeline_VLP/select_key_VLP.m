function clusters = select_key_VLP(clusters,distance_between_VLP,nb_clusters_per_cell, varargin)
    
% selective_kill_parallel_multi_version;
% start_check_parallel_multi_version;

n_clusters  = length(clusters);
x           = zeros(n_clusters,1);
y           = zeros(n_clusters,1);

for i    = 1: n_clusters

   x(i,1)   = mean(clusters(i).xyt(:,1));
   y(i,1)   = mean(clusters(i).xyt(:,2));

end

d      = give_distance_matrix(x,y);
d(d<distance_between_VLP) = 0;

clusters_out(1) = clusters(1);
II              = [1];
indice          = 2;

for i = 2 : min(nb_clusters_per_cell, n_clusters )
    if (sum(d(i,II) < distance_between_VLP ) )==0
       clusters_out(indice) = clusters(i);
       indice = indice + 1;
       II = [II,i];
    end
end

clusters = clusters_out;
clear clusters_out;
    
    
end