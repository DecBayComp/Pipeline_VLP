function [Maps, Maps_densities] = Mapping_without_tracking(movie_per_frame,state,modal,Previous_Assignment,  initial_mesh, name_pipeline, varargin)
warning('off','all');


%% check if multiple processors can be used
% selective_kill_parallel_multi_version;
% start_check_parallel_multi_version;

%% check that data are in microns
[movie_per_frame] = preprocessing(movie_per_frame);

%% ways to use the non tracking
if ~exist('Previous_Assignment'); Previous_Assignment = []; end
if ~exist('initial_mesh');        initial_mesh        = []; end
if ~exist('name_pipeline');       name_pipeline       = 'main'; end

%%
% try
    [~, tout,dt,parameters ] = Assignment_Multi_Mode_Full_Movie(movie_per_frame,modal, Previous_Assignment, name_pipeline);
    [Maps, parameters]       = generate_mesh(tout, parameters, state, initial_mesh);
    [Maps]                   = optimize_Mapping_wihtout_tracking(parameters, Maps, dt, tout,state);
    [Maps]                   = embed_value_meighbors(Maps, state,parameters);
    [Maps]                   = post_processing_local_Mapping_without_tracking(Maps, parameters);
    [Maps]                   = smooth_maps_multi_methods( Maps, 'diffusion' , 'neighbors');
    [Maps_densities]         = output_density_from_Maps(Maps, parameters.desnity_sigma,parameters.density_dx);
    
    %     [Maps, ~, ~, ~]          = generate_density_from_Maps(Maps, parameters.desnity_sigma,parameters.density_dx);
% catch
%     Maps                     = struct;    
% end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  