% RFX stats computed by flipping the sign of subjects R-maps randomly.
%
% Method published in Su et al. (2012) Int. Workshop on Pattern Recognition in NeuroImaging.  
%
% Original author: Li Su 2012-02
% Updated: Cai Wingfield 2015-11, 2016-03
function [observed_map_paths, vertex_level_thresholds] = rfx_tfce(map_paths, n_flips, stat, fdr_thresholds, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.stat.*
    import rsa.util.*
    
    maps_dir = fullfile(userOptions.rootPath, 'Maps');

    n_subjects = numel(userOptions.subjectNames);
    
    
    %% Get actual data
    
    % Load the first dataset to look at size of data
    for chi = 'LR'
        hemi_mesh_stc.(chi) = mne_read_stc_file1(map_paths(1).(chi));
        [n_verts.(chi), n_timepoints] = size(hemi_mesh_stc.(chi).data);
        % delete data fields from hemi_mesh_stc to avoid broadcasting it to all
        % workers
        hemi_mesh_stc.(chi) = rmfield(hemi_mesh_stc.(chi), 'data');
    end
    n_verts_overall = n_verts.L + n_verts.R;
    
    % Compute an adjacency matrix of the downsampled mesh.
    vertex_adjacency = calculateMeshAdjacency(userOptions.targetResolution, userOptions.minDist, userOptions);
    for chi = 'LR'
        % 'iwm' - index within mask
        adjacency_matrix_iwm.(chi) = neighbours2adjacency(hemi_mesh_stc.(chi).vertices, vertex_adjacency);
    end
    
    % Load in and stack up subject correlation-maps
    all_subject_rhos = nan(n_subjects, n_verts_overall, n_timepoints);
    for subject_i = 1:n_subjects
        % Left
        hemi_mesh_stc.L = mne_read_stc_file1(map_paths(subject_i).L);
        all_subject_rhos(subject_i, 1:n_verts.L,       :) = hemi_mesh_stc.L.data;
        % Right
        hemi_mesh_stc.R = mne_read_stc_file1(map_paths(subject_i).R);
        all_subject_rhos(subject_i,   n_verts.L+1:end, :) = hemi_mesh_stc.R.data;
    end
    
    
    %% Simulation    
    
    % We will compute the maximum t-value for each
    % permutation and store those in a null distribution.
    
    % preallocate null distribution vectors for each hemisphere
    h0_l = nan(n_flips, 1);
    h0_r = nan(n_flips, 1);
    
    parfor flip_i = 1:n_flips
        
        % Occasional update
        if mod(flip_i, floor(n_flips/10)) == 0, prints('Flipping coin %d of %d...', flip_i, n_flips); end%if
        
        % Flip a coin for each subject
        flips = (2 * coinToss([n_subjects, 1, 1])) - 1;
        % Copy this to make it the same size as the data
        flips = repmat(flips, [1, n_verts_overall, n_timepoints]);
        
        % Apply the flips to the subject data.
        flipped_rhos = all_subject_rhos .* flips;
        
        if strcmpi(stat, 't')
            % Compute t-stats for this flip
            [h,p,ci, flipped_stats] = ttest(flipped_rhos);

            group_map_sim_both_hemis = squeeze(flipped_stats.tstat);
        elseif strcmpi(stat, 'r')
            % Compute average r-map
            group_map_sim_both_hemis = mean(flipped_rhos, 1);
            group_map_sim_both_hemis = squeeze(group_map_sim_both_hemis);
        else
            error('Must be ''t'' or ''r''.');
        end

        group_map_sim_L = group_map_sim_both_hemis(1:n_verts.L,       :);
        group_map_sim_R = group_map_sim_both_hemis(  n_verts.L+1:end, :);
        
        %% Apply TFCE
        
        group_map_sim_tfce_L = tfce(group_map_sim_L, adjacency_matrix_iwm.L);
        group_map_sim_tfce_R = tfce(group_map_sim_R, adjacency_matrix_iwm.R);
        
        h0_l(flip_i) = max(group_map_sim_tfce_L(:));
        h0_r(flip_i) = max(group_map_sim_tfce_R(:));
    end
    
    h0.L = sort(h0_l);
    h0.R = sort(h0_r);
    
    
    %% Observed maps

    if strcmpi(stat, 't')
        [h,p,ci,stats] = ttest(all_subject_rhos);
        group_map_observed_overall = squeeze(stats.tstat);
    elseif strcmpi(stat, 'r')
        group_map_observed_overall = mean(all_subject_rhos, 1);
        group_map_observed_overall = squeeze(group_map_observed_overall);
        else
            error('Must be ''t'' or ''r''.');
    end
    
    % Split into hemispheres
    group_maps_observed.L = group_map_observed_overall(1:n_verts.L,       :);
    group_maps_observed.R = group_map_observed_overall(  n_verts.L+1:end, :);
    
    for chi = 'LR'
        
        % write out unthresholded map
        observed_map_paths.(chi) = fullfile( ...
            maps_dir, ...
            sprintf('%s_group_%s-map_observed-%sh.stc', userOptions.analysisName, stat, lower(chi)));
        write_stc_file( ...
            hemi_mesh_stc.(chi), ...
            group_maps_observed.(chi), ...
            observed_map_paths.(chi));
        
        % Apply tfce
        group_maps_observed_tfce.(chi) = tfce( ...
            group_maps_observed.(chi), ...
            adjacency_matrix_iwm.(chi));
        
        % write out tfce map
        tfce_map_paths.(chi) = fullfile( ...
            maps_dir, ...
            sprintf('%s_group_%s-map_tfce-%sh.stc', userOptions.analysisName, stat, lower(chi)));
        write_stc_file( ...
            hemi_mesh_stc.(chi), ...
            group_maps_observed_tfce.(chi), ...
            tfce_map_paths.(chi));
        
        % Threshold at corrected p-level
        vertex_level_thresholds.(chi) = nan(size(fdr_thresholds));
        for i = 1:numel(fdr_thresholds)
            vertex_level_thresholds.(chi)(i) = quantile(h0.(chi), 1-fdr_thresholds(i));
        end
    
    end

end%function
