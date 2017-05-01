function [R_mean] = get_average_point_position_Maps(Maps)

n_Maps = length(Maps);
x      = [];
y      = [];

for i = 1: n_Maps;
   x = [x; Maps(i).x];
   y = [y; Maps(i).y];
   
end


if isfield(Maps, 'z');
    z = [];
    for i = 1: n_Maps;
        z = [z; Maps(i).z];
    end

    R_mean(1,1) = mean(x);
    R_mean(2,1) = mean(y);
    R_mean(3,1) = mean(z);
    
    
else
    
    
    R_mean(1,1) = mean(x);
    R_mean(2,1) = mean(y);
    
end






end