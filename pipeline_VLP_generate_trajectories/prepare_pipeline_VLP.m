function [path,files,density_criterion,tt, size_pixel,dt_frame_ms, dt_frame,radius_VLP, state, modal,criterion ,zero_num, nb_clusters_per_cell,  percentage_clusters_per_cell  ,distance_between_VLP,name_pipeline ] = prepare_pipeline_VLP()
%% all parameters for VLP
path                         = pwd;
files                        = subdir(fullfile(pwd, '*.tif'));
t_shift                      = 5;
tt                           = [0:80:1600 - 80];
% tt_shift                   = [t_shift:t_shift:20*t_shift];
density_criterion            = 1000; 
% sigma                      = 0.005;
nb_clusters_per_cell         = 30;
percentage_clusters_per_cell = 0.1;
distance_between_VLP         = 0.05;

size_pixel                   = 0.16;
dt_frame_ms                  = 20;
dt_frame                     = dt_frame_ms/1000;
radius_VLP                   = 1;
state                        = 'diff_voronoi';
modal                        = 'hungarian';
%criterion                    = 'number_point_sliding';
criterion                   = 'fixed_duration_sliding';

zero_num                     = num2str(0);
name_pipeline                = 'vlp';


end