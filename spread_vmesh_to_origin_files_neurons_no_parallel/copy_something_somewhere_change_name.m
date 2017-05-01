function copy_something_somewhere_change_name(original_path, final_path, original_name, final_name)



if ismac||isunix

    original_full = [original_path '/' original_name ];
    final_full    = [final_path    '/' final_name    ];

    
    
%     fprintf('original full %s\n', original_full );
%     fprintf('final    full %s\n', final_full    );
%     
    
    copyfile(original_full, final_full);
    

elseif ispc
    
    original_full = [original_path '\' original_name ];
    final_full    = [final_path    '\' final_name    ];


    copyfile(original_full, final_full);    
    
    
end



end