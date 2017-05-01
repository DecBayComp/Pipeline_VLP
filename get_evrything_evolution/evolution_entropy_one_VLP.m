function [t, S_t]  = evolution_entropy_one_VLP( tau, n_resample)

files   = dir('Maps_*');
n_files = length(files);

t   = nan(n_files,3);
S_t = nan(n_files,1);


for i = 1 : n_files
    inum = num2str(i);
    try
        load(['Maps_' inum '.mat']);
        V               = get_features_out_maps(Maps, 'V');
        V               = V - min(V);
        p               = exp(-V)./sum(exp(-V));
        entropy         = -nansum(p.*log(p));
        t(i,1)          = Maps(1).t_mean;
        t(i,2)          = Maps(1).t_init;
        t(i,3)          = Maps(1).t_end;
        S_t(i,1)        = entropy;
            
    end
    clear Maps;
    
    
end

II  = ~isnan(S_t);
t   = t(II,:);
S_t = S_t(II); 

[t(:,1), S_t] = Gaussian_average_heterogeneous_temporal_data_1D(t(:,1), S_t, tau, n_resample);

Evolution_Entropy.t = t;
Evolution_Entropy.S = S_t;

save('Evolution_Entropy.mat', 'Evolution_Entropy');






end