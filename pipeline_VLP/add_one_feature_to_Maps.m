function [Maps] = add_one_feature_to_Maps(Maps, feature, name_feature)

    n_Maps       = length(Maps);
    
    for i = 1 : n_Maps
        Maps(i).(name_feature) =feature;
    end
    
    

end