% RFX stats computed by flipping the sign of subjects R-maps randomly.
%
% Method published in Su et al. (2012) Int. Workshop on Pattern Recognition in NeuroImaging.  
%
% Original author: Li Su 2012-02
% Updated: Cai Wingfield 2015-11, 2016-03
function [null_t_dists, corrected_thresholds] = cluster_rfx(map_paths, n_flips, p_threshold, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.stat.*
    import rsa.util.*
    
    maps_dir = fullfile(userOptions.rootPath, 'Maps');

    n_subjects = numel(userOptions.subjectNames);
    
    confidence_level = 1 - p_threshold;
    
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
    
    % Write out unthresholded t-map
    for chi = 'LR'
        write_stc_file( ...
            hemi_mesh_stc.(chi), ...
            group_tmaps_observed.(chi), ...
            fullfile(maps_dir, sprintf('%s_group_tmap_observed_%sh.stc', userOptions.analysisName, lower(chi))));
    end
    
    % We will compute the maximum t-value for each
    % permutation and store those in a null distribution.
    
    max_t_values_L = nan(1, n_flips);
    max_t_values_R = nan(1, n_flips);
    max_t_values_B = nan(1, n_flips);
    parfor flip_i = 1:n_flips
        
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

        group_tmap_sim_l = group_tmap_sim(      1:v_l, :);
        group_tmap_sim_r = group_tmap_sim(v_l + 1:end, :);

        max_t_values_L(flip_i) = max(group_tmap_sim_l(:));
        max_t_values_R(flip_i) = max(group_tmap_sim_r(:));
        max_t_values_B(flip_i) = max(group_tmap_sim(:));
    end

    null_t_dists.L = sort(max_t_values_L);
    null_t_dists.R = sort(max_t_values_R);
    null_t_dists.B = sort(max_t_values_B);

    corrected_thresholds.L = quantile(null_t_dists.L, confidence_level);
    corrected_thresholds.R = quantile(null_t_dists.L, confidence_level);
    corrected_thresholds.B = quantile(null_t_dists.L, confidence_level);

end%function

