function [Maps, parameters] = generate_mesh(tout, parameters, state, initial_mesh, varargin)
%%
%

switch  lower(state)

    case {'diff_tree','diff_force_tree', 'diff_drift_tree'}

         [Maps] = build_tree(tout, parameters, initial_mesh);

    case {'diff_voronoi','diff_force_voronoi', 'diff_drift_voronoi'}

         [Maps] = voronoi_mesh( tout, parameters,initial_mesh);

end




end