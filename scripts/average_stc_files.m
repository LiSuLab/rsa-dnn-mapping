function average_stc_paths = average_stc_files(map_paths, name_prefix, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.util.*
    
    mapsDir = fullfile(userOptions.rootPath, 'Maps');
    
    n_files = numel(map_paths);
    
    for chi = 'LR'
        
        for file_i = 1:n_files
            
            prints('Reading stc file %s...', map_paths(file_i).(chi));
            
            stc_struct = mne_read_stc_file1(map_paths(file_i).(chi));
            
            if file_i == 1
                sum_data = stc_struct.data;
            else
                sum_data = sum_data + stc_struct.data;
            end
            
        end%for:file_i
        
        average_data = sum_data / n_files;
        
        % Reuse an stc struct with all the right metadata in it already.
        stc_struct.data = average_data;
        
        average_stc_paths.(chi) = fullfile( ...
            mapsDir, ...
            sprintf('%s-%sh.stc', ...
                name_prefix, ...
                lower(chi)));
        
        prints('Writing averaged data file to %s...', average_stc_paths.(chi));
        
        mne_write_stc_file1(average_stc_paths.(chi), stc_struct);
        
    end%for:chi

end%function
