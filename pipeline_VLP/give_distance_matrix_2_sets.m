function [d] = give_distance_matrix_2_sets(xi, yi,xj, yj,zi,zj, varargin)
%% give the distance matrix 
% to be found in : Matlab/projet/Mapping_without_tracking/prepare_data

if nargin==4
    
    d1 = bsxfun(@plus, xi.^2 , xj'.^2);
    d1 = bsxfun(@minus, d1,   2.* bsxfun(@times, xi, xj') );
    
    d2 = bsxfun(@plus, yi.^2 , yj'.^2);
    d2 = bsxfun(@minus, d2,   2.* bsxfun(@times, yi, yj') );
    
    d = d1 + d2 ;
    
elseif nargin == 6
    
    d1 = bsxfun(@plus, xi.^2 , xj'.^2);
    d1 = bsxfun(@minus, d1,   2.* bsxfun(@times, xi, xj') );
    
    d2 = bsxfun(@plus, yi.^2 , yj'.^2);
    d2 = bsxfun(@minus, d2,   2.* bsxfun(@times, yi, yj') );
    
    d3 = bsxfun(@plus, zi.^2 , zj'.^2);
    d3 = bsxfun(@minus, d3,   2.* bsxfun(@times, zi, zj') );
    
    d = d1 + d2 + d3;
    
end


    
    
    
end