function state = check_state_in_time_evolution(state)


if strcmp(state, 'diff')
   state = 'diff_voronoi'; 
end

if strcmp(state, 'diff_force')
   state = 'diff_force_voronoi'; 
end

if strcmp(state, 'diff_drift')
   state = 'diff_drift_voronoi'; 
end

end