function [test] = check_if_trajectories_already_exist_VLP()

    
try
    test = load('trajectories.txt');
catch
    test = [];
end

end