% RFX stats computed by flipping the sign of subjects R-maps randomly.
%
% Method published in Su et al. (2012) Int. Workshop on Pattern Recognition in NeuroImaging.  
%
% Original author: Li Su 2012-02
% Updated: Cai Wingfield 2015-11, 2016-03
function [observed_map_paths, permuted_map_paths] = rfx_permutation_tmaps(map_paths, n_flips, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.stat.*
    import rsa.util.*
    
    maps_dir = fullfile(userOptions.rootPath, 'Maps');

    n_subjects = numel(userOptions.subjectNames);
    
    % Load the first dataset to look at size of data
    for chi = 'LR'
        hemi_mesh_stc.(chi) = mne_read_stc_file1(map_paths(1).(chi));
        [n_verts.(chi), n_timepoints] = size(hemi_mesh_stc.(chi).data);
    end
    n_verts_overall = n_verts.L + n_verts.R;
    
    % Preallocate
    all_subject_rhos = nan(n_subjects, n_verts_overall, n_timepoints);
    
    for subject_i = 1:n_subjects
        % Left
        hemi_mesh_stc.L = mne_read_stc_file1(map_paths(subject_i).L);
        all_subject_rhos(subject_i, 1:n_verts.L,       :) = hemi_mesh_stc.L.data;
        % Right
        hemi_mesh_stc.R = mne_read_stc_file1(map_paths(subject_i).R);
        all_subject_rhos(subject_i,   n_verts.L+1:end, :) = hemi_mesh_stc.R.data;
    end

    [h,p,ci,stats] = ttest(all_subject_rhos);
    group_tmap_observed_overall = squeeze(stats.tstat);
    
    % Set nan values to 0
    % TODO: why would there be nans?
    group_tmap_observed_overall(isnan(group_tmap_observed_overall)) = 0;
    
    % Split into hemispheres
    group_tmaps_observed.L = group_tmap_observed_overall(1:n_verts.L,       :);
    group_tmaps_observed.R = group_tmap_observed_overall(  n_verts.L+1:end, :);
    
    % Write out unthresholded t-maps
    for chi = 'LR'
        observed_map_paths.(chi) = fullfile( ...
            maps_dir, ...
            sprintf('%s_group_tmap_observed-%sh.stc', userOptions.analysisName, lower(chi)));
        write_stc_file( ...
            hemi_mesh_stc.(chi), ...
            group_tmaps_observed.(chi), ...
            observed_map_paths.(chi));
    end
    
    % We will compute the maximum t-value for each
    % permutation and store those in a null distribution.
    
    % -1 because we'll 0-index them
    n_digits = ceil(log(n_flips-1));
    flip_format = ['%0', num2str(n_digits), 'd'];
    
    % delete data fields from hemi_mesh_stc to avoid broadcasting it to all
    % workers
    for chi = 'LR'
        hemi_mesh_stc.(chi) = rmfield(hemi_mesh_stc.(chi), 'data');
    end
    
    parfor flip_i = 1:n_flips
        
        flip_name = sprintf(flip_format, flip_i-1);
        
        % Occasional update
        if mod(flip_i, floor(n_flips/100)) == 0, prints('Flipping coin %d of %d...', flip_i, n_flips); end%if
        
        % Flip a coin for each subject
        flips = coinToss([n_subjects, 1, 1]);
        % Copy this to make it the same size as the data
        flips = repmat(flips, [1, n_verts_overall, n_timepoints]);
        
        % Apply the flips to the subject data.
        flipped_rhos = all_subject_rhos .* flips;
        
        % Compute t-stats for this flip
        [h,p,ci, flipped_stats] = ttest(flipped_rhos);
        
        group_tmap_sim = squeeze(flipped_stats.tstat);

        group_tmap_sim_l = group_tmap_sim(1:n_verts.L,       :);
        group_tmap_sim_r = group_tmap_sim(  n_verts.L+1:end, :);
        
        % Write out permuted tmap
        permuted_map_paths(flip_i).L = fullfile( ...
            maps_dir, ...
            sprintf('rfx_perm%s_tmap-lh.stc', flip_name));
        permuted_map_paths(flip_i).R = fullfile( ...
            maps_dir, ...
            sprintf('rfx_perm%s_tmap-rh.stc', flip_name));
        
        write_stc_file( ...
            hemi_mesh_stc.L, ...
            group_tmap_sim_l, ...
            permuted_map_paths(flip_i).L);
        write_stc_file( ...
            hemi_mesh_stc.R, ...
            group_tmap_sim_r, ...
            permuted_map_paths(flip_i).R);
    end

end%function

