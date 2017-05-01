function start_check_parallel_multi_version

status_all = 'all';


Version = version ;
if ismac||isunix
    numero  = str2num(Version(13:16));
elseif ispc
    numero  = str2num(Version(16:19));
end

if (numero <= 2013 )
    if ( matlabpool('size') == 0 )
        myCluster = parcluster('local');
        switch status_all
            case 'half' 
                nb_workers = floor(myCluster.NumWorkers./2);
            case 'all'
                nb_workers = (myCluster.NumWorkers);
        end
        matlabpool(myCluster, nb_workers);
        pctRunOnAll  warning('off','all'); 
    end
else
    if( isempty(gcp('nocreate')) )
        myCluster = parcluster('local');
        switch status_all
            case 'all' 
                nb_workers = (myCluster.NumWorkers);
            case 'half'
                nb_workers = floor(myCluster.NumWorkers./2);
        end
        parpool(nb_workers );
        pctRunOnAll  warning('off','all'); 
    end
end