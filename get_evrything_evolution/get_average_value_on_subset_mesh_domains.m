function average_value = get_average_value_on_subset_mesh_domains(Maps, indice_in, features)


n_in          = length(indice_in);
value         = [];

for i = 1 : n_in
    
    value = [ value; Maps( indice_in(i) ).(features) ];
    
end


average_value = nanmean(value);



end