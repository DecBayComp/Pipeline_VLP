function Maps = optimize_Mapping_wihtout_tracking(parameters, Maps, dt, tout,state, varargin)
%% optimize various fields in mapping withtout tracking


switch lower(state)

    %%
    case 'diff_tree'

        [Maps] = loop_D(parameters, tout, Maps,dt);
    %%    
    case 'diff_force_tree'

        [Maps] = loop_D_f(parameters, tout, Maps,dt);
    %%    
    case 'diff_drift_tree'

        [Maps] = loop_D_v(parameters, tout, Maps,dt);
    %%                
    case 'diff_voronoi'

        [Maps] = give_D(parameters,Maps,dt);
    %%   
    case 'diff_force_voronoi'

        [Maps] = give_D_f(parameters,  Maps,dt);
    %%   
    case 'diff_drift_voronoi'

        [Maps] = give_D_v(parameters,  Maps,dt);


end    



end