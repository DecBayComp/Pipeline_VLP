function get_evrything_evolution(target_folder, varargin)
warning off;
%%
% if ischar(i)
%      i = str2num(i);
% end

if nargin>= 1
    cd(target_folder);
end

original_path    = pwd;
extension        = 'vmesh';
[liste, n_liste] = search_files_hierarchy_individual_folders(target_folder,extension) ;

%% parameters done for VLP
radius     = 0.05;
tau        = 120;
n_resample = 1e4;

%%
for i = 1 : n_liste
    fprintf('%i\n', i);
    get_evolution_diffusivity(liste, n_liste,i, target_folder,radius,tau, n_resample);
    get_evolution_diffusivity_variability(liste, n_liste,i,target_folder,radius,tau, n_resample);
    get_evolution_entropy(liste, n_liste,i, target_folder, tau, n_resample);
    get_evolution_free_energy(liste, n_liste,i, target_folder, tau, n_resample);
    get_evolution_potential_depth(liste, n_liste,i, target_folder,  tau, n_resample);
end
%%

cd (original_path);

end