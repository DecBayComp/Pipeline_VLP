function parameters = generate_parameters(movie_per_frame, parameters,name_pipeline,  varargin)
%% generate all parameters for analysis
% to be found in : Matlab/projet/Mapping_without_tracking/parameters

if ~exist('name_pipeline');        name_pipeline        = 'main'; end

switch lower(name_pipeline)
    
    case 'main'

        parameters = generate_parameters_main(movie_per_frame, parameters);
    
    case 'vlp'
        
        parameters = generate_parameters_vlp(movie_per_frame, parameters);
        
    case 'rac1'
        
        parameters = generate_parameters_rac1(movie_per_frame, parameters);
        
    case 'nucleus'
        
        parameters = generate_parameters_nucleus(movie_per_frame, parameters);
        
    case 'neurons'
        
        parameters = generate_parameters_neurons(movie_per_frame, parameters);
        
    case 'large_maps'
        
        parameters = generate_parameters_large_maps(movie_per_frame, parameters);
        
    case 'gaba'
        
        parameters = generate_parameters_gaba(movie_per_frame, parameters);
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
















