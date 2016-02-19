% Assumes that fiels are in increasing order of lags.
%
% Cai Wingfield 2016-02
function average_stc_paths = average_stc_files(map_paths, name_prefix, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.util.*
    
    mapsDir = fullfile(userOptions.rootPath, 'Maps');
    
    n_lags = numel(map_paths);
    
    for chi = 'LR'
        
        for lag_i = 1:n_lags
            
            prints('Reading stc file %s...', map_paths(lag_i).(chi));
            
            stc_struct = mne_read_stc_file1(map_paths(lag_i).(chi));
            
            this_tmin = stc_struct.tmin;
            this_tmax = stc_struct.tmax;
            
            if lag_i == 1
                sum_data = stc_struct.data;
        
                sum_data_tmin = stc_struct.tmin;
                sum_data_tmax = stc_struct.tmax;
                % This should always be the same
                tstep = stc_struct.tstep;
            else
                
                % WE ARE ASSUMING THAT EACH INCOMING FILE IS EXACTLY tstep
                % BEHIND THE PREVIOUS.
                
                % We need to truncate the beginning of this data
                tdiff = sum_data_tmin - this_tmin;
                % debug
                assert(tdiff == tstep);
                
                % Trim the first frame of this data
                this_data = stc_struct.data(2:end, :);
                
                % Trim the last frame of the data pool
                sum_data = sum_data(1:end-1, :);
                sum_data_tmax = this_tmax;
                
                sum_data = sum_data + this_data;
            end
            
        end%for:file_i
        
        % Mean the remaining data
        average_data = sum_data / n_lags;
        
        % Reuse an stc struct with all the right vertices and tstep in it already.
        stc_struct.data = average_data;
        
        % Correct the the limits
        stc_struct.tmin = sum_data_tmin;
        stc_struct.tmax = sum_data_tmax;
        
        % Write it out
        average_stc_paths.(chi) = fullfile( ...
            mapsDir, ...
            sprintf('%s-%sh.stc', name_prefix, lower(chi)));
        
        prints('Writing averaged data file to %s...', average_stc_paths.(chi));
        
        mne_write_stc_file1(average_stc_paths.(chi), stc_struct);
        
    end%for:chi

end%function
