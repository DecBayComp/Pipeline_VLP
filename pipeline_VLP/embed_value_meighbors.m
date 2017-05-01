function [Maps] = embed_value_meighbors(Maps,state,parameters)

switch lower(state)
    
    case {'diff_tree', 'diff_force_tree', 'diff_drift_tree'}
    
        if parameters.d == 2
            
            Maps = embed_value_meighbors_tree_2D(Maps);
            
        elseif parameters.d == 3
            
            Maps = embed_value_meighbors_tree_3D(Maps);
        end
        
        
    case {'diff_voronoi','diff_force_voronoi','diff_drift_voronoi'}
        
        if parameters.d == 2
            
            Maps = embed_value_meighbors_voronoi_2D(Maps);
            
        elseif parameters.d == 3

            Maps = embed_value_meighbors_voronoi_3D(Maps);
            
        end
        
end

end