function threshold = rfx_threshold(rfx_paths, confidence_level)

    import rsa.*
    import rsa.meg.*
    
    % convert a 0.05 into a 0.95;
    quantile_level = 1 - confidence_level;
    
    n_perms = numel(rfx_paths);
    
    for chi = 'LR'
        
        null_pool = [];
        for file_i = 1:n_perms
        
            prints('\tReading stc file for perm %d...', rfx_paths(file_i).(chi));
            stc_struct = mne_read_stc_file1(rfx_paths(file_i).(chi));
            
            % Add the values for this perm to the null pool
            null_pool = [null_pool; stc_struct.data(:)];
        end
        
        null_pool = sort(null_pool);
        
        threshold.(chi) = quantile(null_pool, quantile_level);
        
    end%for:chi

end
