function movie_per_frame = check_units_data(movie_per_frame)




%% number of dimentions

if isfield(movie_per_frame, 'z')
     d                         = 3;
else
    d                         = 2;
end

    
xx = [];
yy = [];

if d == 2
    for i = 1 : length(movie_per_frame)    
        xx = [xx; movie_per_frame(i).x];
        yy = [yy; movie_per_frame(i).y];
    end
elseif d ==3
    zz = [];
    for i = 1 : length(movie_per_frame)    
        xx = [xx; movie_per_frame(i).x];
        yy = [yy; movie_per_frame(i).y];
        zz = [zz; movie_per_frame(i).z]; 
    end
end

xx_min = min(xx);
xx_max = max(xx);

yy_min = min(yy);
yy_max = max(yy);

if d==3
    zz_min = min(zz);
    zz_max = max(zz);
end

if d==2
    l_tot = sqrt(  (xx_max - xx_min).^2 + (yy_max - yy_min).^2);
elseif d==3
    l_tot = sqrt(  (xx_max - xx_min).^2 + (yy_max - yy_min).^2 + (zz_max - zz_min).^2);
end


if (l_tot < 1e-3)
   
    if d ==2
        for i = 1 : length(movie_per_frame)   
            movie_per_frame(i).x = movie_per_frame(i).x*1e6;
            movie_per_frame(i).y = movie_per_frame(i).y*1e6;
        end
    elseif d==3
        for i = 1 : length(movie_per_frame)   
            movie_per_frame(i).x = movie_per_frame(i).x*1e6;
            movie_per_frame(i).y = movie_per_frame(i).y*1e6;
            movie_per_frame(i).z = movie_per_frame(i).z*1e6;
        end
    end
    
    
end



end