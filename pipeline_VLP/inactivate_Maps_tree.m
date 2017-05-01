function Maps = inactivate_Maps_tree(Maps, parameters)


D_init = parameters.D_init;
for i = 1 : length(Maps)
   if ( Maps(i).D == D_init)
       Maps(i).D = NaN;
   end
end


end