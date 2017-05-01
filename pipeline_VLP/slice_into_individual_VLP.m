function [clusters]  = slice_into_individual_VLP(xyt, x_target, y_target, radius )

% selective_kill_parallel_multi_version;
% start_check_parallel_multi_version;

n_x_target = length(x_target);
for i = 1 : n_x_target
    xx = x_target(i); yy = y_target(i);
%     r2              = (xyt(:,1) - xx).^2 + (xyt(:,2) - yy).^2 ;
%     JJ              = r2 <= radius.^2;
    JJ =  ( xyt(:,1) >= xx - radius  ) & ( xyt(:,1) <= xx + radius) ...
        & ( xyt(:,2) >= yy - radius  ) & ( xyt(:,2) <= yy + radius);

    clusters(i).xyt = xyt(JJ,:);
    clusters(i).nb  = length(clusters(i).xyt(:,1));
end

nb_tot = [];
for  i = 1 : n_x_target
    nb_tot = [nb_tot ; clusters(i).nb];
end

[~,I] = sort(nb_tot,'descend');

cluster2 = clusters;
for i =1 : n_x_target
    cluster2(i).xyt = clusters(I(i)).xyt ;
    cluster2(i).nb  = clusters(I(i)).nb  ;
end

clusters = cluster2;
clear cluster2;



if ~exist('clusters')
    clusters = [];
end

end