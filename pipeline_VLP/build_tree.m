function [mesh]  = build_tree(tout, parameters, mesh_init, varargin)
%% Generate the mesh quadtree in 2D octree in 3D
%to be found in : Matlab/projet/Mapping_without_tracking/tree


if parameters.d == 2
    [mesh]  = build_tree_2D(tout, parameters, mesh_init);
elseif parameters.d ==3
    [mesh]  = build_tree_3D(tout, parameters, mesh_init);
end
 

end