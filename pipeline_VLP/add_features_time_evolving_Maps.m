function Maps = add_features_time_evolving_Maps(Maps,t_init, t_end )



n = length(Maps);
for i = 1 : n
    
   Maps(i).t_init = t_init;
   Maps(i).t_end  = t_end;
   Maps(i).t_mean = mean([t_init t_end]);
    
end






end