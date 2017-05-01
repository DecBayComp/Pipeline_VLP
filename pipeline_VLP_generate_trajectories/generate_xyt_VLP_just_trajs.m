function [xyt, x_target, y_target, t_min, t_max] = generate_xyt_VLP_just_trajs(test,name_folder,filename, pathname,size_pixel, dt_frame_ms,density_criterion,zero_num,name_pipeline )



if isempty(test)
    xyt                  = video_to_localisations_VLP(filename, pathname,size_pixel, dt_frame_ms);
%     xyt                  = shift_time_VLP(xyt,tt , tt_shift);
    xyt(:,3)             = regularize_time_VLP(xyt(:,3));
    t_min                = min(xyt(:,3));
    t_max                = max(xyt(:,3));
%     [x_target, y_target] = get_target_densities_VLP(xyt, density_criterion);
    movie_per_frame      = convert_traj_struct(xyt); 
    cd ..
    trajs                = regenerate_and_clean_trajectories_VLP(movie_per_frame , zero_num,name_pipeline);
    cd(name_folder);
    clear xyt
    xyt(:,1) = trajs(:,2); xyt(:,2) = trajs(:,3); xyt(:,3) = trajs(:,4);
    x_target = []; y_target = [];
else
    trajs = test;
    xyt(:,1) = trajs(:,2); xyt(:,2) = trajs(:,3); xyt(:,3) = trajs(:,4);
    t_min                = min(xyt(:,3));
    t_max                = max(xyt(:,3));
    %[x_target, y_target] = get_target_densities_VLP(xyt, density_criterion);
    x_target = []; y_target= [];
end






end