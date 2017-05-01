function parameters = get_param_BP_parameters(parameters)


%% 
parameters.kind_init                     = 'random';
parameters.all_methods_init              = {'zeros','random', 'normal', 'mix'};


%% damping factor for update of the messages
parameters.lambda                        = 0.75;

%% loops for updates
parameters.number_stop                   = 10;
parameters.global_loop                   = 1;
parameters.nb_duration_state_1           = 100;  
parameters.number_stop_loop_tree         = 5;



end