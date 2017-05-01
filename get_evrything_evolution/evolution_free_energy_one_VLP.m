function [t, F_t]  = evolution_free_energy_one_VLP(tau, n_resample)

files   = dir('Maps_*');
n_files = length(files);

t   = nan(n_files,3);
F_t = nan(n_files,1);


for i = 1 : n_files
    inum = num2str(i);
    try
        load(['Maps_' inum '.mat']);
        V               = get_features_out_maps(Maps, 'V');
        V               = V - min(V);
        F               = -log( nansum(exp(-V)) )- ( -log(length(V)));
        t(i,1)          = Maps(1).t_mean;
        t(i,2)          = Maps(1).t_init;
        t(i,3)          = Maps(1).t_end;
        F_t(i,1)        = F;
            
    end
    clear Maps;
    
    
end
II  = ~isnan(F_t);
t   = t(II,:);
F_t = F_t(II); 

[t(:,1), F_t] = Gaussian_average_heterogeneous_temporal_data_1D(t(:,1), F_t, tau, n_resample);


Evolution_Free_Energy.t = t;
Evolution_Free_Energy.F = F_t;

save('Evolution_Free_Energy.mat', 'Evolution_Free_Energy');






end