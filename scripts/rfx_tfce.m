% RFX stats computed by flipping the sign of subjects R-maps randomly.
%
% Method published in Su et al. (2012) Int. Workshop on Pattern Recognition in NeuroImaging.  
%
% Original author: Li Su 2012-02
% Updated: Cai Wingfield 2015-11, 2016-03
function [observed_map_paths] = rfx_tfce(map_paths, n_flips, stat, fdr_threshold, userOptions)

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
        if mod(flip_i, floor(n_flips/100)) == 0, prints('Flipping coin %d of %d...', flip_i, n_flips); end%if
        
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
        
        group_map_sim_L = tfce(group_map_sim_L, adjacency_matrix_iwm.L);
        group_map_sim_R = tfce(group_map_sim_R, adjacency_matrix_iwm.R);
        
        h0_l(flip_i) = max(group_map_sim_L);
        h0_r(flip_i) = max(group_map_sim_R);
    end
    
    h0.L = sort(h0_l); clear h0_l;
    h0.R = sort(h0_r); clear h0_r;
    
    
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
        group_maps_observed_tfce.(chi) = tfce(group_maps_observed.(chi));
        
        % Threshold at corrected p-level
        vertex_level_threshold = quantile(h0.(chi), 1-fdr_threshold);
        
        % Threshold to reveal supra-corrected-p-threshold vertices
        corrected_thresholded_maps.(chi) = zeros(size(group_maps_observed_tfce.(chi)));
        corrected_thresholded_maps.(chi)(corrected_thresholded_maps.(chi) > vertex_level_threshold) = corrected_thresholded_maps.(chi)(corrected_thresholded_maps.(chi) > vertex_level_threshold);
        
        % write out corrected cluster map
        observed_thresholded_map_paths_corr.(chi) = fullfile( ...
            maps_dir, ...
            sprintf('%s_group_%s-map_observed_thresholded_corrected-%sh.stc', userOptions.analysisName, stat, lower(chi)));
        write_stc_file( ...
            hemi_mesh_stc.(chi), ...
            corrected_thresholded_maps.(chi), ...
            observed_thresholded_map_paths_corr.(chi));
    
    end

end%function


%% %%%%%%%%%%%%%
% Subfunctions %
%%%%%%%%%%%%%%%%

% vertex_adjacency - we assume that this has already been downsampled to the current resolution.
% vertices - a list of vertex names (as opposed to vertex indices) relating
%            to the current data.
function adjacency_matrix_iwm = neighbours2adjacency(masked_vertices, vertex_adjacency)

    import rsa.*
    import rsa.meg.*
    
    n_vertices = numel(masked_vertices);
    
    % To produce an adjacency matrix is to enumerate all the edges of a
    % graph. So we will describe each edge by a "starting" and an "ending"
    % vertex.
    
    % For each starting vertex, there are at most `max_valence` corresponding ending
    % vertices, where `max_valence` is the maximum valence of the graph.
    max_valence = size(vertex_adjacency, 2);
    
    % Before we proceed, we will forget all adjacency information about
    % vertices outside the mask.
    % 'iwm' - index within mask
    vertex_adjacency_iwm = vertex_adjacency(masked_vertices, :);
    
    % The row indices of vertex_adjacency_iwm now point to vertices in
    % masked_vertices.  Hence we call it iwm.  However the actual ENTRIES
    % in vertex_adjacency_iwm are vertices indexed-within-brain.
    
    % The list of starting vertices is the list of vertices repeated
    % `max_valence` times.
    % 'iwb' - index within brain
    starting_vertices_iwb = repmat(masked_vertices(:), max_valence, 1); 
    
    % Since vertex_adjacency is indexed in the first dimension by actual
    % vertex name, the corresponding list of ending vertices can be
    % achieved by concatenating the rows of vertex_adjacency.
    ending_vertices_iwb = vertex_adjacency_iwm;
    ending_vertices_iwb = ending_vertices_iwb(:);
    
    % Since not every vertex achieves `max_valence`, we look for nans in
    % ending vertices and then delete these from both lists
    null_edge_locations = isnan(ending_vertices_iwb);
    starting_vertices_iwb = starting_vertices_iwb(~null_edge_locations);
    ending_vertices_iwb   = ending_vertices_iwb(~null_edge_locations);
    
    % We finally need to check again for ending vertices outside the mask,
    % and remove those.
    inside_mask_edge_locations = ismember(ending_vertices_iwb, masked_vertices);
    starting_vertices_iwb = starting_vertices_iwb(inside_mask_edge_locations);
    ending_vertices_iwb   = ending_vertices_iwb(inside_mask_edge_locations);
    
    % Our adjacency matrix will descibe vertices inside the mask only, so
    % the row and columns of the matrix will be pointers into the list of
    % vertices in the mask.
    % 'iwm' - index within mask
    starting_vertices_iwm = nan(1, numel(starting_vertices_iwb));
    ending_vertices_iwm   = nan(1, numel(ending_vertices_iwb));
    for vi = 1:numel(starting_vertices_iwb)
        % lookup within-mask indices of within-brain vertex names
        starting_vertices_iwm(vi) = find(masked_vertices == starting_vertices_iwb(vi), 1);
        ending_vertices_iwm(vi) = find(masked_vertices == ending_vertices_iwb(vi), 1);
    end
    
    % Now we create the adjacency matrix from the edge lists.
    % Need to convert to doubles for no good reason.
    adjacency_matrix_iwm = sparse(double(starting_vertices_iwm), double(ending_vertices_iwm), 1, n_vertices, n_vertices);

end

function cluster_stats = calculate_cluster_stats(labelled_spatiotemporal_clusters, group_maps, vertex_level_threshold)

...%cluster_stats = cluster_peak( ...
...%cluster_stats = cluster_extent( ...
cluster_stats = cluster_exceedence_mass( ...
            labelled_spatiotemporal_clusters, ...
            group_maps, ...
            vertex_level_threshold);

end
