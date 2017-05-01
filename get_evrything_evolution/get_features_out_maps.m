function vector = get_features_out_maps(Maps, features)

n_Maps = length(Maps);
vector = [];

for i = 1 : n_Maps;
    
    vector = [vector; Maps(i).(features)];
    
end






end