function [t, D_t]  = evolution_diffusivity_one_VLP(radius, tau, n_resample)

files   = dir('Maps_*');
n_files = length(files);

t   = nan(n_files,3);
D_t = nan(n_files,1);


for i = 1 : n_files
    inum = num2str(i);
    try
        load(['Maps_' inum '.mat']);
        [R_mean]      = get_average_point_position_Maps(Maps);
        indice_in     = get_mesh_index_inside_sphere_region(Maps, R_mean, radius);
        D             = get_average_value_on_subset_mesh_domains(Maps, indice_in, 'D');
        t(i,1)        = Maps(1).t_mean;
        t(i,2)        = Maps(1).t_init;
        t(i,3)        = Maps(1).t_end;
        D_t(i,1)      = D;
        
    end
    clear Maps;
    
    
end

II  = ~isnan(D_t);
t   = t(II,:);
D_t = D_t(II); 

[t(:,1), D_t] = Gaussian_average_heterogeneous_temporal_data_1D(t(:,1), D_t, tau, n_resample);

Evolution_Diffusivite.t = t;
Evolution_Diffusivite.D = D_t;

save('Evolution_Diffusivite.mat', 'Evolution_Diffusivite');






end