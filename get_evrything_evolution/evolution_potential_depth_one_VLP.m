function [t, depth_t]  = evolution_potential_depth_one_VLP(tau, n_resample)

files   = dir('Maps_*');
n_files = length(files);

t   = nan(n_files,3);
depth_t = nan(n_files,1);


for i = 1 : n_files
    inum = num2str(i);
    try
        load(['Maps_' inum '.mat']);

        V             = get_features_out_maps(Maps, 'V');
        indice_out    = get_out_domain_mesh_Maps(Maps);
        min_V         = min(V);
        max_V         = mean(V(indice_out));
        depth         = max_V - min_V;
        t(i,1)        = Maps(1).t_mean;
        t(i,2)        = Maps(1).t_init;
        t(i,3)        = Maps(1).t_end;
        depth_t(i,1)  = depth;
        
    end
    clear Maps;
    
    
    
end

II  = ~isnan(depth_t);
t   = t(II,:);
depth_t =depth_t(II); 

[t(:,1), depth_t] = Gaussian_average_heterogeneous_temporal_data_1D(t(:,1),depth_t, tau, n_resample);

Evolution_Depth.t     = t;
Evolution_Depth.depth = depth_t;

save('Evolution_Depth.mat', 'Evolution_Depth');






end