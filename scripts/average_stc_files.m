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
            
            if lag_i == 1
                sum_data = stc_struct.data;
                reusable_stc_struct = stc_struct;
                
                % To record the number of model rs that have contributed to
                % this datapoint, for the purposes of correctly averaging
                % later.
                % ones() because we've got at least one datapoint from this
                % zero-lag model.
                n_timepoints = size(sum_data, 2);
                n_contributions = ones(1, n_timepoints);
            else
                
                % Trim the first frame of this data
                this_data = stc_struct.data(2:end, :);
                
                % Record this contribution
                width_this_data = size(this_data, 1);
                n_contributions(1:width_this_data) = n_contributions(1:width_this_data) + 1;
                
                sum_data = sum_data + this_data;
            end
            
        end%for:file_i
        
        % Mean the remaining data, appropriate
        average_data = sum_data ./ repmat(n_contributions, n_verts, 1);
        
        % Reuse an stc struct with all the right vertices and tstep in it already.
        reusable_stc_struct.data = average_data;
        
        % Write it out
        average_stc_paths.(chi) = fullfile( ...
            mapsDir, ...
            sprintf('%s-%sh.stc', name_prefix, lower(chi)));
        
        prints('Writing averaged data file to %s...', average_stc_paths.(chi));
        
        mne_write_stc_file1(average_stc_paths.(chi), reusable_stc_struct);
        
    end%for:chi

end%function
