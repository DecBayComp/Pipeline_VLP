function Maps = switch_status_state(Maps, new_status)
%% switch state status base on local size and diffusion
% Matlab/projet/Mapping_without_tracking/tree


for i = 1 : length(Maps)
    Maps(i).status = new_status;
end
    
end