function [x_start, y_start, densities, XX, YY, radius] = set_seed_densities_multi_VLP(trajs, value)


x                           = trajs(:,2);
y                           = trajs(:,3);

[densities, XX, YY]         = effective_densities(x, y, 0.1);
[x_start, y_start, II, JJ]  = give_center_VLP(trajs, value);
n_start                     = length(x_start);

if isempty(x_start)
    radius = [];
else
    

    for i = 1 : n_start
%         fprintf('%i\n', i);
        max_density       = densities(JJ(i),II(i));
        max_density_4     = max_density./4;
        low               = 0.1*max_density;

        densities_loc     = densities;
        KK                = densities_loc(:) <= low;
        densities_loc(KK) = 0;

        x                 = XX(:);
        y                 = YY(:);
%         epsilon           = 0.01*max_density;
        epsilon           = 0.02*max_density;

        d                 = sqrt( (x - x_start(i)).^2 + (y - y_start(i)).^2) ;
        [d,LL]            = sort(d,'ascend');
        densities_ordered = densities_loc(LL);
        
        
        MM               = find(  (densities_ordered > max_density_4 - epsilon/2 ) & ( densities_ordered < max_density_4 + epsilon/2 ) ); 
        radius(i,1)       = d(MM(1));
        clear LL d  KK densities_loc;

    end

end



end