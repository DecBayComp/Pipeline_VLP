function d = give_distance_matrix(x,y,z, varargin)
%% give the distance matrix from a set of points
%to be found in : Matlab/projet/Mapping_without_tracking/prepare_data


if nargin ==2
    
    d1 = bsxfun(@plus, x.^2 , x'.^2);
    d1 = bsxfun(@minus, d1,   2.* bsxfun(@times, x, x') );
    
    d2 = bsxfun(@plus, y.^2 , y'.^2);
    d2 = bsxfun(@minus, d2,   2.* bsxfun(@times, y, y') );
    
    d = d1 + d2;
    
elseif nargin ==3
    
    d1 = bsxfun(@plus, x.^2 , x'.^2);
    d1 = bsxfun(@minus, d1,   2.* bsxfun(@times, x, x') );
    
    d2 = bsxfun(@plus, y.^2 , y'.^2);
    d2 = bsxfun(@minus, d2,   2.* bsxfun(@times, y, y') );
    
    d3 = bsxfun(@plus, z.^2 , z'.^2);
    d3 = bsxfun(@minus, d3,   2.* bsxfun(@times, z, z') );
    
    
    d = d1 + d2 + d3;
    
end

   



    