% Assumes that fiels are in increasing order of lags.
%
% Cai Wingfield 2016-02
function average_stc_paths = lag_integrate_stc_files(map_paths, name_prefix, userOptions, truncate)

    import rsa.*
    import rsa.meg.*
    import rsa.util.*
    
    mapsDir = fullfile(userOptions.rootPath, 'Maps');
    
    n_lags = numel(map_paths);
    
    if ~exist('truncate', 'var'), truncate = false; end
    
    for chi = 'LR'
        
        for lag_i = 1:n_lags
            
            prints('(%d/%d) Reading stc file %s...', lag_i, n_lags, map_paths(lag_i).(chi));
            
            stc_struct = mne_read_stc_file1(map_paths(lag_i).(chi));
            
            if lag_i == 1
                reusable_stc_struct = stc_struct;
                
                % To record the number of model rs that have contributed to
                % this datapoint, for the purposes of correctly averaging
                % later.
                % ones() because we've got at least one datapoint from this
                % zero-lag model.
                [n_verts, n_timepoints] = size(stc_struct.data);
                
                % Preallocate
                aggregate_data = zeros(n_verts, n_timepoints, n_lags);
                contribution_counts = zeros(1, n_timepoints);
            end
            
            trimmed_epoch_length = n_timepoints - lag_i + 1;
            
            aggregate_data(:, 1:trimmed_epoch_length, lag_i) = stc_struct.data(:, lag_i:end);
            
            contribution_counts(1:trimmed_epoch_length) = contribution_counts(1:trimmed_epoch_length) + 1;
            
        end%for:lag_i
        
        % Mean the remaining data, appropriate
        average_data = sum(aggregate_data, 3) ./ repmat(contribution_counts, n_verts, 1);
        
        % Truncate those parts of the data after which there is sub-optimal
        % SNR
        if truncate
            trimmed_epoch_length = n_timepoints - n_lags + 1;
            average_data = average_data(:, 1:trimmed_epoch_length);
        end
        
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
