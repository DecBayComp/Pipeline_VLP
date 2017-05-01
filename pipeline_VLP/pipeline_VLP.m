function [ clusters, xyt, x_target, y_target,size_pixel,dt_frame] = pipeline_VLP(i, target_folder_start_process, varargin)
warning off;   

if ischar(i)
     i = str2num(i);
end
cd(target_folder_start_process);

%% to avoid weird stuff
destroy_files_evrywhere_hierarchy(target_folder_start_process,'*.full');
destroy_files_evrywhere_hierarchy(target_folder_start_process,'*.2d');
%%
[path,files,density_criterion,tt, size_pixel,dt_frame_ms,...
    dt_frame,radius_VLP, state, modal,criterion ,zero_num ,...
   nb_clusters_per_cell,  percentage_clusters_per_cell, distance_between_VLP, ...
   name_pipeline  ] = prepare_pipeline_VLP();

n_files = length(files);
%%
if i <= n_files

    % for i = 1 : length(files)
%%
    [name_folder, filename,pathname, inum] = dynamical_folder_and_path_generation_VLP(i, files);
    [test] = check_if_trajectories_already_exist_VLP();
    cd(name_folder);
%%
    [xyt, x_target, y_target, t_min, t_max] = generate_xyt_VLP(test,name_folder, ...
        filename, pathname,size_pixel, dt_frame_ms,density_criterion,zero_num ,name_pipeline) ;
%%
    if ~isempty(x_target)
        [clusters]           = slice_into_individual_VLP(xyt, x_target, y_target, radius_VLP ); 
        [clusters]           = select_key_VLP(clusters, distance_between_VLP, nb_clusters_per_cell);
        [clusters]           = cluster_to_movie_per_frame_VLP(clusters, t_min, t_max, dt_frame);
    %plot_clusters_VLP( clusters);
%%
        for j = 1 : length(clusters)   
            fprintf('%i\n', j);
            jnum            = num2str(j); 
            movie_per_frame = clusters(j).movie_per_frame;
            trajs           = regenerate_trajectories_VLP(movie_per_frame , jnum, name_pipeline);
            time_evolving_maps_short(jnum,movie_per_frame,state,modal, criterion,name_pipeline);
            generate_cluster_files_VLP(jnum);
        end
    end
        
end
        
% end

cd(path);
% collect_all_clusters_files_VLP()



% close all force;
% clear;clc;
if  ~exist('clusters')
   clusters = []; 
end
if  ~exist('xyt')
   xyt = []; 
end

%selective_kill_parallel_multi_version;


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

