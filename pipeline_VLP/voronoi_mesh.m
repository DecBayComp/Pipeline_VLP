function mesh = voronoi_mesh( tout, parameters,initial_mesh, varargin)
%% voronoi mesh generation
%to be found in : Matlab/projet/Mapping_without_tracking/voronoi

if parameters.d == 2
    [ mesh ] = voronoi_meshing( tout, parameters,initial_mesh );
elseif parameters.d ==3 
    [ mesh ] = voronoi_meshing_3D( tout, parameters,initial_mesh );
end




end