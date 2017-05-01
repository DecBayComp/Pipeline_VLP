function [ clusters, xyt, x_target, y_target,size_pixel,dt_frame] = pipeline_VLP_generate_trajectories(i, target_folder, varargin)
warning off;   

if ischar(i)
    i = str2num(i);
end
cd(target_folder);

%% to avoid weird stuff
destroy_files_evrywhere_hierarchy(target_folder, '*.full');
destroy_files_evrywhere_hierarchy(target_folder, '*.2d');

[path,files,density_criterion,tt, size_pixel,dt_frame_ms,...
    dt_frame,radius_VLP, state, modal,criterion ,zero_num ,...
   nb_clusters_per_cell,  percentage_clusters_per_cell, distance_between_VLP,name_pipeline ] = prepare_pipeline_VLP();

n_files = length(files);

if i<= n_files

% for i = 1 : length(files)
        
    [name_folder, filename,pathname, inum] = dynamical_folder_and_path_generation_VLP(i, files);
    [test] = check_if_trajectories_already_exist_VLP();
    cd(name_folder);
    
    [xyt, x_target, y_target, t_min, t_max] = generate_xyt_VLP_just_trajs(test,name_folder, ...
        filename, pathname,size_pixel, dt_frame_ms,density_criterion,zero_num,name_pipeline ) ;
    
end

        
% end

% cd(path);

% close all force;
% clear;clc;
if  ~exist('clusters')
   clusters = []; 
end
if  ~exist('xyt')
   xyt = []; 
end

% selective_kill_parallel_multi_version;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

