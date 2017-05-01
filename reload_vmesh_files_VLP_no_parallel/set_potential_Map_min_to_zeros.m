function Maps = set_potential_Map_min_to_zeros(Maps)


nn = length(Maps);
V  = zeros(nn,1);
for i = 1 : nn
    V(i) = Maps(i).V;
end
V = V - min(V);

for i = 1 : nn
    Maps(i).V = V(i);
end

end