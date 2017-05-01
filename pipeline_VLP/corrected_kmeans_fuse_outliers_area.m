function [IDX, C, parameters] = corrected_kmeans_fuse_outliers_area(R, initial_C, parameters)


number_of_cluster   = parameters.number_of_cluster;
min_number_per_zone = parameters.min_number_per_zone;

[IDX, C, ~, D]  = kmeans(R, number_of_cluster, 'Start',initial_C, 'EmptyAction', 'singleton','OnlinePhase','off','Options',parameters.opts);
% compteur        = 0;
% compteur_lim    = number_of_cluster;


% fprintf('%i\n', number_of_cluster);
indice = 0;
while(1)
    [Id,  IDX, C,D ]       = detect_undersampled_mesh_domains_and_remove_one(IDX,C,D, min_number_per_zone );
    if isempty(Id)
        break;
    end
    indice = indice + 1;
%     fprintf('%i\n', indice);
end

[n_C, m_C]                   = size(C);
parameters.number_of_cluster = n_C;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [Id,  IDX, C,D ] = detect_undersampled_mesh_domains_and_remove_one(IDX,C,D, min_number_per_zone)

[n,m]      = size(IDX);
[n_C, m_C] = size(C);
[n_D, m_D] = size(D);
nb         = zeros( n_C,1);
II         = [1 : n_C]';

for i      = 1 : n_C
    JJ     = IDX == i;
    nb(i)  = sum(JJ);
end

KK        = nb < min_number_per_zone;
Id        = II(KK);
clear     KK ;
n_Id      = length(Id);

clear II JJ KK;
if isempty(Id)
else
%     i           = 1;
    i           = randi(n_Id,1,1);
    JJ          = Id(i);
    II          = IDX == JJ;
    
    D(II,JJ)    = nan;
    [Y, KK]     = min(D, [], 2);
    IDX(II)     = KK(II);
    
    III         = unique(IDX);
    I4          = [1:length(III)];
    IDX_out     = zeros(n,1);
    C_out       = zeros(n_C-1,m_C);
    D_out       = zeros(n,n_C-1);
    
    for j = 1 : length(III)
       J_loc           = III(j);
%        fprintf('j %i\t J_loc %i\n', j, J_loc);
       II              = IDX == J_loc;
       IDX_out(II)     = I4(j);
       C_out(j,:)      = C(J_loc,:); 
       D_out(:,j)      = D(:,J_loc) ;
    end
    C                 = C_out   ;
    IDX               = IDX_out ; 
    D                 = D_out   ;
end
% unique(IDX)

% number_of_cluster = number_of_cluster-1;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Id = detect_undersampled_mesh_domains(IDX, number_of_cluster,min_number_per_zone )

[n,m] = size(IDX);
nb    = zeros( number_of_cluster,1);
II    = [1:number_of_cluster]';
for i = 1 : number_of_cluster
    JJ = IDX == i;
    nb(i) = sum(JJ);
end

KK = nb < min_number_per_zone;
Id = II(KK);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IDX, C] = reassign(R, Id, IDX, number_of_cluster);



for i = 1 : length(Id)
   JJ          = Id(i);
   II          = IDX == JJ;
   D_loc       = D(II,:);
   D_loc(:,JJ) = nan; 
   [Y, KK]     = min(D_loc, [], 2);
   IDX(II)     = KK; 
   
    
end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% while (compteur <= compteur_lim)
%     Id = detect_undersampled_mesh_domains(IDX, number_of_cluster,min_number_per_zone );
%     if isempty(Id)
%         break;
%     else
%         [n,m]             = size(Id);
%         II                = randi(n,1,1);
%         Id_kill           = Id(II);
%         C(Id_kill,:)      = [];
%         number_of_cluster = number_of_cluster - 1;
%         
%         
%         
% %         [IDX, C, ~]       = kmeans(R, number_of_cluster, 'Start',C, 'EmptyAction', 'singleton','OnlinePhase','off','Options', parameters.opts);
% %         [initial_C]       = initialize_voronoi_random(R, number_of_cluster);
% %         fprintf('number of cluster %i\n',number_of_cluster );
% %         [initial_C]       = initialize_voronoi_bubble(R, number_of_cluster); 
% %         [nn,mm]           = size(initial_C);
% %         fprintf('number_cluster %i\t column %i\t rows %i\n', number_of_cluster , mm, nn);
% %         [IDX, C, ~]       = kmeans(R, number_of_cluster, 'Start',initial_C, 'EmptyAction', 'singleton','OnlinePhase','off','Options', parameters.opts);
%         clear initial_C;
% 
%     end
%     
%     
%     compteur = compteur + 1; 
% %     fprintf('%i\n', compteur);
% end