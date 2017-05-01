function get_evolution_entropy(liste, n_liste,i, target_folder,tau, n_resample, varargin)


% if nargin>= 1
%    cd(target_folder); 
% end
% radius           = 0.05;
original_path    = pwd;
% extension        = 'vmesh';
% [liste, n_liste] = search_files_hierarchy_individual_folders(target_folder,extension) ;

% for i = 1 : n_liste
     if (i<= n_liste)  
   try 
        cd(liste{i});
        evolution_entropy_one_VLP(tau, n_resample);
        
   end
     end
% end


cd(original_path);



end