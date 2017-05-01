function [indice_range, n_window] = get_indice_range(criterion,movie_per_frame, name_pipeline, varargin)

if ~exist('name_pipeline');       name_pipeline       = 'main'; end

parameters          = generate_parameters(movie_per_frame,[], name_pipeline);
duration            = parameters.duration   ;
t_sliding           = parameters.t_sliding  ;
nb_points           = parameters.nb_points  ;
nb_sliding          = parameters.nb_sliding ;
t_init              = parameters.t_init ;      
t_end               = parameters.t_end;      
n_movie_per_frame   = parameters.n_movie_per_frame ;


criterion = lower(criterion);

switch criterion
    %%
    case 'fixed_duration_sequence'

        n_window = floor((t_end - t_init)./duration);
        t_tot    = zeros(n_movie_per_frame, 1);
        for i = 1 : n_movie_per_frame
            t_tot(i,1) = movie_per_frame(i).t;
        end
        t_tot = t_tot - t_tot(1);

        indice_range = zeros(n_window,2);
        for i = 1 : n_window            
            II                = find( (t_tot <= i*duration) & (t_tot >= (i-1)*duration) );       
            indice_range(i,1) = II(1);
            indice_range(i,2) = II(end);   
        end
        
        if n_window == 0
            indice_range(1,1) = 1;
            indice_range(1,2) = n_movie_per_frame;
        end

    case 'fixed_duration_sliding'
        
        n_window = floor((t_end - t_init - duration)./t_sliding + 1);
        t_tot    = zeros(n_movie_per_frame, 1);
        for i = 1 : n_movie_per_frame
            t_tot(i,1) = movie_per_frame(i).t;
        end
        t_tot = t_tot - t_tot(1);
        
        indice_range = zeros(n_window,2);
        for i = 1 : n_window 
            II                = find( (t_tot <= (i-1)*t_sliding + duration) & (t_tot >= (i-1)*t_sliding) );       
            indice_range(i,1) = II(1);
            indice_range(i,2) = II(end);
        end
        
        
        if n_window == 0
            indice_range(1,1) = 1;
            indice_range(1,2) = n_movie_per_frame;
        end
        
    case 'number_point_sequence'
        
        
        nb_tot    = zeros(n_movie_per_frame-1, 1);
        for i = 1 : n_movie_per_frame-1
            nb_tot(i,1) = movie_per_frame(i).nb_linkable;
        end
        nb_tot_cumul    = zeros(n_movie_per_frame-1, 1);
        nb_tot_cumul(1,1) = nb_tot(1,1);
        for i = 2 :n_movie_per_frame-1
            nb_tot_cumul(i,1) = nb_tot_cumul(i-1,1) + nb_tot(i,1);
        end
        
        n_window = floor(nb_tot_cumul(end)./nb_points);
        indice_range = zeros(n_window,2);
        for i = 1 : n_window 
            II                = find( (nb_tot_cumul <= i*nb_points) & (nb_tot_cumul >= (i-1)*nb_points) );       
            indice_range(i,1) = II(1);
            indice_range(i,2) = II(end);
        end
                
        if n_window == 0
            indice_range(1,1) = 1;
            indice_range(1,2) = n_movie_per_frame-1;
        end
        
    case 'number_point_sliding'
        
        nb_tot    = zeros(n_movie_per_frame-1, 1);
        for i = 1 : n_movie_per_frame-1
            nb_tot(i,1) = movie_per_frame(i).nb_linkable;
        end
        nb_tot_cumul    = zeros(n_movie_per_frame-1, 1);
        nb_tot_cumul(1,1) = nb_tot(1,1);
        for i = 2 :n_movie_per_frame-1
            nb_tot_cumul(i,1) = nb_tot_cumul(i-1,1) + nb_tot(i,1);
        end
        n_window     = floor( (nb_tot_cumul(end,1) - nb_points)./nb_sliding + 1);
        indice_range = zeros(n_window,2);
        for i = 1 : n_window  
            II                = find( (nb_tot_cumul <= (i-1)*nb_sliding + nb_points) & (nb_tot_cumul >= (i-1)*nb_sliding) );       
            indice_range(i,1) = II(1);
            indice_range(i,2) = II(end);
        end
        
        
        if n_window == 0
            indice_range(1,1) = 1;
            indice_range(1,2) = n_movie_per_frame-1;
        end
end