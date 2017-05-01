function std_value = get_std_value_on_subset_mesh_domains(Maps, indice_in, features)


n_in          = length(indice_in);
value         = [];

for i = 1 : n_in
    
    value = [ value; Maps( indice_in(i) ).(features) ];
    
end


std_value = nanstd(value);



end