function Maps = Set_potential_stochastic_integral(Maps, lambda, sigma, dt_theo, varargin)

if isfield(Maps,'V')

    n_Maps = length(Maps);
    V      = zeros(n_Maps,1);
    D      = zeros(n_Maps,1);

    for i = 1 : n_Maps
        V(i,1) = Maps(i).V;
        D(i,1) = Maps(i).D;
    end

    II = D<0.;
    JJ = isnan(D);

    if  (sum(II)==0) && (sum(JJ)==0)

        V = V + lambda*log(D);


    elseif (sum(II)>0) && (sum(JJ)==0)

        D(II) = sigma*sigma/dt_theo;
        V     = V + lambda*log(D);

    else


    end

    for i = 1 : n_Maps
        Maps(i).V = V(i,1);
        Maps(i).D = D(i,1) ;
    end

    Maps = set_potential_Map_min_to_zeros(Maps);
else
    
end

end