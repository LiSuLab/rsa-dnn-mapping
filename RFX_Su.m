% RFX stats computed by flipping the sign of subjects R-maps randomly.
%
% Method published in Su et al. (2012) Int. Workshop on Pattern Recognition in NeuroImaging.  
%
% Original author: Li Su 2012-02
% Updated: Cai Wingfield 2015-11
function vertex_level_threshold = RFX_Su(map_paths, STC_metadatas, n_flips, confidence_level, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.util.*

    n_subjects = numel(userOptions.subjectNames);
    
    % Load the first dataset to look at size of data
    [v_l, n_timepoints] = size(directLoad(map_paths(1).L));
    [v_r, n_timepoints] = size(directLoad(map_paths(1).R));
    n_verts = v_l + v_r;
    
    prints('RFX ...');
    
    all_rhos = nan(n_subjects, n_verts, n_timepoints);
    
    for subject_i = 1:n_subjects
        % Left
        all_rhos(subject_i,       1:v_l, :) = directLoad(map_paths(subject_i).L);
        % Right
        all_rhos(subject_i, v_l + 1:end, :) = directLoad(map_paths(subject_i).R);
    end

    [h,p,ci,stats] = ttest(all_rhos);
    true_group_ts = squeeze(stats.tstat);
    
    % Set nan values to 0
    true_group_ts(isnan(true_group_ts)) = 0;
    
    % Split into hemispheres
    true_group_ts_l = true_group_ts(      1:v_l, :);
    true_group_ts_r = true_group_ts(v_l + 1:end, :);
    
    % Write out unthresholded t-maps
    unthresholded_t_path = fullfile(userOptions.rootPath, 'Meshes', 't_unthresholded');
    write_stc_file(STC_metadatas(1).L, true_group_ts_l, [unthresholded_t_path, '-lh.stc']);
    write_stc_file(STC_metadatas(1).R, true_group_ts_r, [unthresholded_t_path, '-rh.stc']);
    
    max_t_values = zeros(1, n_flips);

    for flip_i = 1:n_flips
        if mod(flip_i, floor(n_flips/100)) == 0, prints('Flipping coin %d of %d...', flip_i, n_flips); end%if
        
        % Flip a coin for each subject
        flips = ((rand(n_subjects, 1, 1) > 0.5) .* 2) - 1;
        % Copy this to make it the same size as the data
        flips = repmat(flips, [1, n_verts, n_timepoints]);
        
        % Apply the flips to the subject data.
        flipped_rhos = all_rhos .* flips;
        
        % Compute t-stats for this flip
        [h,p,ci, flipped_stats] = ttest(flipped_rhos);
        
        flipped_ts = squeeze(flipped_stats.tstat);

        flipped_ts_l = flipped_ts(      1:v_l, :);
        flipped_ts_r = flipped_ts(v_l + 1:end, :);

        max_t_values(flip_i) = max(flipped_ts(:));
    end

    null_t_distribution = sort(max_t_values);

    vertex_level_threshold = null_t_distribution(ceil(size(null_t_distribution, 2) * confidence_level));
    
    % Threshold data
    
    % Split into hemispheres
    thresholded_group_ts_l = true_group_ts_l;
    thresholded_group_ts_l(thresholded_group_ts_l < vertex_level_threshold) = 0;
    thresholded_group_ts_r = true_group_ts_r;
    thresholded_group_ts_r(thresholded_group_ts_r < vertex_level_threshold) = 0;
    
    % Write out thresholded t-maps
    thresholded_t_path = fullfile(userOptions.rootPath, 'Meshes', 't_thresholded');
    write_stc_file(STC_metadatas(1).L, thresholded_group_ts_l, [thresholded_t_path, '-lh.stc']);
    write_stc_file(STC_metadatas(1).R, thresholded_group_ts_r, [thresholded_t_path, '-rh.stc']);

end%function

