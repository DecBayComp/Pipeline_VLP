function [Maps] = give_D(parameters,Maps,dt, varargin)

% selective_kill_parallel_multi_version;
% start_check_parallel_multi_version;


n_zone                          = length(Maps);
sigma                           = parameters.sigma;
opt                             = parameters.opt;
D_out                           = zeros(n_zone,1);
Feval                           = zeros(n_zone,1);
if parameters.d == 2

    for i = 1 : n_zone  
     %    fprintf('optimizing D in %i\n', i);
        dr2    = Maps(i).dx.^2 + Maps(i).dy.^2;
        if isempty(dr2)
            D_out(i) = NaN;
            Feval(i) = -inf;
        else
            D_init = sum(dr2/(4.*dt))./length(dr2);
            Dmax    = 2*D_init;
            Dmin    = D_init./5;
            [D_out(i), Feval(i)]                 = fminbnd(@(D) likelihood_D(dr2,dt, D,sigma),Dmin, Dmax, opt); 
        end
    end


    for i = 1 : n_zone;

        Maps(i).D                = D_out(i);
        Maps(i).log_likelihood   = -Feval(i);
    end

elseif parameters.d ==3

    for i = 1 : n_zone  
        % fprintf('optimizing D in %i\n', i);
        dr2    = Maps(i).dx.^2 + Maps(i).dy.^2 + Maps(i).dz.^2;
        if isempty(dr2)
            D_out(i) = NaN;
            Feval(i) = -inf;
        else
            D_init = sum(dr2/(6.*dt))./length(dr2);
            Dmax    = 2*D_init;
            Dmin    = D_init./5;
            [D_out(i), Feval(i)]                 = fminbnd(@(D) likelihood_D_3D(dr2,dt, D,sigma),Dmin, Dmax, opt); 
        end
    end


    for i = 1 : n_zone;

        Maps(i).D                = D_out(i);
        Maps(i).log_likelihood   = -Feval(i);
    end


end
        
        
        
        
        
        
        

end