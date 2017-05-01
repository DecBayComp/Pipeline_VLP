function kill_parallel_multi_version



Version = version ;
if ismac||isunix
    numero  = str2num(Version(13:16));
elseif ispc
    numero  = str2num(Version(16:19));
end

if (numero <= 2013 )
    if ( matlabpool('size') ~= 0 )
        matlabpool close force local;
    end
else
    if( ~isempty(gcp('nocreate')) )
        delete(gcp('nocreate'));
    end
end







