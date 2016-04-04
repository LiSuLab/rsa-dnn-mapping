% RFX stats computed by flipping the sign of subjects R-maps randomly.
%
% Method published in Su et al. (2012) Int. Workshop on Pattern Recognition in NeuroImaging.  
%
% Original author: Li Su 2012-02
% Updated: Cai Wingfield 2015-11, 2016-03
function [observed_map_paths, corrected_ps] = rfx_cluster(map_paths, n_flips, primary_p_threshold, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.stat.*
    import rsa.util.*
    
    maps_dir = fullfile(userOptions.rootPath, 'Maps');
    simulation_dir = fullfile(userOptions.rootPath, 'Sim');

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
        adjacency_matrix_iwm.(chi) = neighbours2adjacency(hemi_mesh_stc.(chi).vertices, vertex_adjacency, userOptions);
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
        
        % Compute t-stats for this flip
        [h,p,ci, flipped_stats] = ttest(flipped_rhos);
        
        group_tmap_sim_both_hemis = squeeze(flipped_stats.tstat);

        group_tmap_sim_L = group_tmap_sim_both_hemis(1:n_verts.L,       :);
        group_tmap_sim_R = group_tmap_sim_both_hemis(  n_verts.L+1:end, :);
        
        % For some reason Matlab won't let me do this loop insid of a
        % parfor.
        
        chi = 'L';
        
        [labelled_sim_clusters, sim_cluster_stats] = compute_cluster_stats(adjacency_matrix_iwm.(chi), group_tmap_sim_L, primary_p_threshold);
        
        h0_l(flip_i) = max(sim_cluster_stats);
        
        chi = 'R';
        
        [labelled_sim_clusters, sim_cluster_stats] = compute_cluster_stats(adjacency_matrix_iwm.(chi), group_tmap_sim_R, primary_p_threshold);
        
        h0_r(flip_i) = max(sim_cluster_stats);
    end
    
    h0.L = h0_l;
    h0.R = h0_r;
    
    clear h0_l h0_r;
    
    
    %% Observed t-maps

    [h,p,ci,stats] = ttest(all_subject_rhos);
    group_tmap_observed_overall = squeeze(stats.tstat);
    
    % Set nan values to 0
    % TODO: why would there be nans?
    group_tmap_observed_overall(isnan(group_tmap_observed_overall)) = 0;
    
    % Split into hemispheres
    group_tmaps_observed.L = group_tmap_observed_overall(1:n_verts.L,       :);
    group_tmaps_observed.R = group_tmap_observed_overall(  n_verts.L+1:end, :);
    
    
    %% Identify observed clusters
    
    for chi = 'LR'
        
        [labelled_spatiotemporal_clusters.(chi), cluster_stats.(chi)] = compute_cluster_stats(adjacency_matrix_iwm.(chi), group_tmaps_observed.(chi), primary_p_threshold);
        
        % write out unthresholded t-map
        observed_map_paths.(chi) = fullfile( ...
            maps_dir, ...
            sprintf('%s_group_tmap_observed-%sh.stc', userOptions.analysisName, lower(chi)));
        write_stc_file( ...
            hemi_mesh_stc.(chi), ...
            group_tmaps_observed.(chi), ...
            observed_map_paths.(chi));
        
        % write out cluster map
        cluster_labels_map_paths.(chi) = fullfile( ...
            m, ...
            sprintf('%s_group_tmap_observed_clusters-%sh.stc', userOptions.analysisName, lower(chi)));
        write_stc_file( ...
            hemi_mesh_stc.(chi), ...
            labelled_spatiotemporal_clusters.(chi), ...
            cluster_labels_map_paths.(chi));
        
    end%for:chi
    
    
    %% Assign corrected ps to each observed cluster
    
    for chi = 'LR'
        
        n_clusters = numel(unique(labelled_spatiotemporal_clusters.(chi)));
        
        for cluster_i = 1:n_clusters
            corrected_ps.(chi)(cluster_i) = quantile(cluster_stats.(chi)(cluster_i), h0.(chi));
        end
        
    end

end%function


%% %%%%%%%%%%%%%
% Subfunctions %
%%%%%%%%%%%%%%%%

% vertex_adjacency - we assume that this has already been downsampled to the current resolution.
% vertices - a list of vertex names (as opposed to vertex indices) relating
%            to the current data.
function adjacency_matrix_iwm = neighbours2adjacency(masked_vertices, vertex_adjacency, userOptions)

    import rsa.*
    import rsa.meg.*
    
    n_vertices = numel(masked_vertices);
    
    % To produce an adjacency matrix is to enumerate all the edges of a
    % graph. So we will describe each edge by a "starting" and an "ending"
    % vertex.
    
    % Before we proceed, we will forget all adjacency information about
    % vertices outside the mask
    vertex_adjacency = vertex_adjacency(masked_vertices, :);
    
    % For each starting vertex, there are at most `max_valence` corresponding ending
    % vertices, where `max_valence` is the maximum valence of the graph.
    max_valence = size(vertex_adjacency, 2);
    
    % The list of starting vertices is the list of vertices repeated
    % `max_valence` times.
    % 'iwb' - index within brain
    starting_vertices_iwb = repmat(masked_vertices(:), max_valence, 1); 
    
    % Since vertex_adjacency is indexed in the first dimension by actual
    % vertex name, the corresponding list of ending vertices can be
    % achieved by concatenating the rows of vertex_adjacency.
    ending_vertices_iwb = vertex_adjacency';
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

% Compute the connected components of a graph based on a sparse adjacency
% matrix.
%
% via http://blogs.mathworks.com/steve/2007/03/20/connected-component-labeling-part-3/
function component_list = connected_components(adjacency_matrix)

    if numel(adjacency_matrix) == 0
       component_list = {}; 
    else
    
        % Add 1s to the diagonal of the adjacency matrix to ensure that each
        % lone vertex becomes a connected component.
        adjacency_matrix(1:size(adjacency_matrix, 1)+1:end) = 1;

        % Compute the Dulmage-Mendelsohn decomposition of the adjacency matrix.
        [row_perm, col_perm, row_blockdiv, col_blockdiv] = dmperm(adjacency_matrix);

        % Then `row_blockdiv` contains the indices in `row_perm` beginning each
        % connected component in the adjacency matrix. -1 for fenceposting.
        n_connected_components = numel(row_blockdiv)-1;

        component_list = cell(n_connected_components, 1);
        for comp_i = 1:n_connected_components
           component_list{comp_i} = row_perm(row_blockdiv(comp_i):row_blockdiv(comp_i+1)-1);
        end
    end
end


function spatial_cluster_labels = label_spatiotemporal_clusters(thresholded_tmap, adjacency_matrix)

    [n_verts, n_timepoints] = size(thresholded_tmap);
    
    % Don't want to remember size info about this.
    adjacency_matrix = full(adjacency_matrix);

	% To begin with we don't look for temporal contiguity.  We label
    % spatially contiguous clusters in each timepoint with a unique
    % label.
    spatial_cluster_labels = zeros(n_verts, n_timepoints);
    
    running_cluster_count = 0;
    for t = 1:n_timepoints
        
        % We're interested in identifying contiguous clusters, so we'll
        % forget adjacency information for sub-threshold vertices.
        
        % We just take the sub-matrix for the vertices above the threshold
        % 'iwm' - index within mask
        super_threshold_vs_iwm = find(thresholded_tmap(:, t));
        % 'iwc' - index within cluster
        thresholded_adjacency_matrix_iwc = adjacency_matrix(super_threshold_vs_iwm, super_threshold_vs_iwm);
        
        % A cell array of components at this timepoint. Each component is a vector of
        % vertex names.
        component_list_this_t_iwc = connected_components(thresholded_adjacency_matrix_iwc);
        n_clusters_this_t = numel(component_list_this_t_iwc);
        
        for component_i = 1:n_clusters_this_t
            running_cluster_count = running_cluster_count + 1;
            
            % Change the vis in the connected components to point back into the
            % mask region.
            % 'iwc' - index within cluster
            component_vs_iwc = component_list_this_t_iwc{component_i};
            component_vs_iwm = super_threshold_vs_iwm(component_vs_iwc);
            
            spatial_cluster_labels(component_vs_iwm, t) = running_cluster_count;
        end
    end
    
    merging_done_last_pass = true;
    while merging_done_last_pass
        % If merging was done last pass, we'll check agian.
        
        merging_done_last_pass = false;
    
        % Make a forward pass and relabel any clusters which are
        % temporally adjacent

        for t = 1:n_timepoints-1
            vis_overlap = find( ...
                spatial_cluster_labels(:, t) .* spatial_cluster_labels(:, t+1));

            % record pairs of cluster ids which are to be merged between
            % these adjacent timepoints
            cluster_id_pairs_to_merge = [];
            for overlap_vi = vis_overlap'
                cluster_pair = [spatial_cluster_labels(overlap_vi, t), spatial_cluster_labels(overlap_vi, t+1)];
                % only remember to merge genuinely different pairs
                if cluster_pair(1) ~= cluster_pair(2)
                    cluster_id_pairs_to_merge = [ ...
                        cluster_id_pairs_to_merge; ...
                        cluster_pair];
                end
            end
            cluster_id_pairs_to_merge = unique(cluster_id_pairs_to_merge, 'rows');
            
            % Remember if merging will be done this pass
            if numel(cluster_id_pairs_to_merge) > 0
                merging_done_last_pass = true;
            end

            % merge clusters
            for cluster_id_pair = cluster_id_pairs_to_merge'

                % relabel all vertices of the matching clusters in both
                % layers to have the same label
                spatial_cluster_labels(spatial_cluster_labels == max(cluster_id_pair)) = min(cluster_id_pair);

                % relabel all to-merge clusters too
                cluster_id_pairs_to_merge(cluster_id_pairs_to_merge == max(cluster_id_pair)) = min(cluster_id_pair);
            end
        end%for

        % Now make a reverse pass

        for t = n_timepoints:-1:2
            vis_overlap = find( ...
                spatial_cluster_labels(:, t) .* spatial_cluster_labels(:, t-1));

            % record pairs of cluster ids which are to be merged between
            % these adjacent timepoints
            cluster_id_pairs_to_merge = [];
            for overlap_vi = vis_overlap'
                cluster_pair = [spatial_cluster_labels(overlap_vi, t), spatial_cluster_labels(overlap_vi, t-1)];
                % only remember to merge genuinely different pairs
                if cluster_pair(1) ~= cluster_pair(2)
                    cluster_id_pairs_to_merge = [ ...
                        cluster_id_pairs_to_merge; ...
                        cluster_pair];
                end
            end
            cluster_id_pairs_to_merge = unique(cluster_id_pairs_to_merge, 'rows');
            
            if numel(cluster_id_pairs_to_merge) > 0
                merging_done_last_pass = true;
            end

            % merge clusters
            for cluster_id_pair = cluster_id_pairs_to_merge'

                % relabel all vertices of the matching clusters in both
                % layers to have the same label
                spatial_cluster_labels(spatial_cluster_labels == max(cluster_id_pair)) = min(cluster_id_pair);

                % relabel all to-merge clusters too
                cluster_id_pairs_to_merge(cluster_id_pairs_to_merge == max(cluster_id_pair)) = min(cluster_id_pair);
            end
        end%for
    end%while
    
    % Relabel clusters
    
    % these will be sorted
    remaining_cluster_labels = unique(spatial_cluster_labels);
    
    % forget the 'background' cluster
    remaining_cluster_labels = remaining_cluster_labels(remaining_cluster_labels > 0);
    
    n_clusters = numel(remaining_cluster_labels);
    
    % renumber clusters in ascending order to prevent collisions
    for rem_cluster_lab_i = 1:n_clusters
        rem_cluster_lab = remaining_cluster_labels(rem_cluster_lab_i);
        spatial_cluster_labels(spatial_cluster_labels == rem_cluster_lab) = rem_cluster_lab_i;
    end
end

function [labelled_spatiotemporal_clusters, cluster_stats] = compute_cluster_stats(adjacency_matrix, group_tmaps_observed, primary_p_threshold)
    vertex_level_threshold = quantile(group_tmaps_observed(:), 1-primary_p_threshold);
    thresholded_tmap = (group_tmaps_observed > vertex_level_threshold);
    labelled_spatiotemporal_clusters = label_spatiotemporal_clusters(thresholded_tmap, adjacency_matrix);
    
    n_clusters = numel(unique(labelled_spatiotemporal_clusters));
    
    cluster_stats = nan(n_clusters, 1);
    for cluster_i = 1:n_clusters
        % cluster exceedence mass
        this_cluster_location = (labelled_spatiotemporal_clusters == cluster_i);
        cluster_exceedences = group_tmaps_observed(this_cluster_location) - vertex_level_threshold;
        cluster_stats(cluster_i) = sum(cluster_exceedences(:));
    end
end

