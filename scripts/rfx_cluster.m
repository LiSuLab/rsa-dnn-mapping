% RFX stats computed by flipping the sign of subjects R-maps randomly.
%
% Method published in Su et al. (2012) Int. Workshop on Pattern Recognition in NeuroImaging.  
%
% Original author: Li Su 2012-02
% Updated: Cai Wingfield 2015-11, 2016-03
function [observed_map_paths, permuted_map_paths] = rfx_cluster(map_paths, n_flips, primary_p_threshold, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.stat.*
    import rsa.util.*
    
    maps_dir = fullfile(userOptions.rootPath, 'Maps');
    simulation_dir = fullfile(userOptions.rootPath, 'Sim');

    n_subjects = numel(userOptions.subjectNames);
    
    % Compute an adjacency matrix of the downsampled mesh.
    vertex_adjacency = calculateMeshAdjacency(userOptions.targetResolution, userOptions.minDist);
    
    % Load the first dataset to look at size of data
    for chi = 'LR'
        hemi_mesh_stc.(chi) = mne_read_stc_file1(map_paths(1).(chi));
        [n_verts.(chi), n_timepoints] = size(hemi_mesh_stc.(chi).data);
    end
    n_verts_overall = n_verts.L + n_verts.R;
    
    
    %% Observed t-maps
    
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
    
    
    %% Simulation
    
    
    
    
    %% Identify observed clusters
    
    for chi = 'LR'
        
        vertex_level_threshold = quantile(group_tmaps_observed.(chi), 1-primary_p_threshold);
        thresholded_tmap = (group_tmaps_observed.(chi) > vertex_level_threshold);
    
        adjacency_matrix = sparse_connectivity_matrix_from_vert_adjacency(vertices, vertex_adjacency);
        
        labelled_spatiotemporal_clusters.(chi) = label_spatiotemporal_clusters(thresholded_tmap, adjacency_matrix);
        n_clusters = numel(unique(labelled_spatiotemporal_clusters.(chi)));
            
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
            simulation_dir, ...
            sprintf('%s_group_tmap_observed_clusters-%sh.stc', userOptions.analysisName, lower(chi)));
        write_stc_file( ...
            hemi_mesh_stc.(chi), ...
            labelled_spatiotemporal_clusters.(chi), ...
            cluster_labels_map_paths.(chi));
        
        for cluster_i = 1:n_clusters
            % cluster exceedence mass
            vertices_this_cluster = (labelled_spatiotemporal_clusters.(chi) == cluster_i);
            cluster_exceedences = group_tmaps_observed.(chi)(vertices_this_cluster, :) - vertex_level_threshold;
            cluster_stats.(chi)(cluster_i) = sum(cluster_exceedences(:));
        end
        
        
    
    end%for:chi
    
    
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
            simulation_dir, ...
            sprintf('rfx_perm%s_tmap-lh.stc', flip_name));
        permuted_map_paths(flip_i).R = fullfile( ...
            simulation_dir, ...
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

% vertex_adjacency - we assume that this has already been downsampled to the current resolution.
% vertices - a list of vertex names (as opposed to vertex indices) relating
%            to the current data.
function adjacency_matrix = sparse_connectivity_matrix_from_vert_adjacency(vertices, vertex_adjacency)

    import rsa.*
    import rsa.meg.*
    
    % To produce an adjacency matrix is to enumerate all the edges of a
    % graph. So we will describe each edge by a "starting" and an "ending"
    % vertex.
    
    % For each starting vertex, there are at most `max_valence` corresponding ending
    % vertices, where `max_valence` is the maximum valence of the graph.
    max_valence = size(vertex_adjacency, 2);
    
    % The list of starting vertices is the list of vertices repeated
    % `max_valence` times.
    starting_vertices = repmat(vertices(:), max_valence, 1); 
    
    % Since vertex_adjacency is indexed in the first dimension by actual
    % vertex name, the corresponding list of ending vertices can be
    % achieved by concatenating the rows of vertex_adjacency.
    ending_vertices = vertex_adjacency';
    ending_vertices = ending_vertices(:);
    
    % Since not every vertex achieves `max_valence`, we look for nans in
    % ending vertices and then delete these from both lists
    null_edge_locations = isnan(ending_vertices);
    starting_vertices = starting_vertices(~null_edge_locations);
    ending_vertices   = ending_vertices(~null_edge_locations);
    
    % Now we create the adjacency matrix from the edge lists.
    adjacency_matrix = sparse(starting_vertices, ending_vertices, 1);

end

% Compute the connected components of a graph based on a sparse adjacency
% matrix.
%
% via http://blogs.mathworks.com/steve/2007/03/20/connected-component-labeling-part-3/
function component_list = connected_components(sp_adjacency_matrix)
    
    % Add 1s to the diagonal of the adjacency matrix to ensure that each
    % lone vertex becomes a connected component.
    sp_adjacency_matrix(1:size(sp_adjacency_matrix, 1):end);
    
    % Compute the Dulmage-Mendelsohn decomposition of the adjacency matrix.
    [row_perm, col_perm, row_blockdiv, col_blockdiv] = dmperm(sp_adjacency_matrix);
    
    % Then `row_blockdiv` contains the indices in `row_perm` beginning each
    % connected component in the adjacency matrix.
    n_connected_components = numel(row_blockdiv)-1;
    
    for comp_i = 1:n_connected_components
       component_list{comp_i} = row_perm(row_blockdiv(comp_i):row_blockdiv(comp_i+1));
    end
end


function spatial_cluster_labels = label_spatiotemporal_clusters(thresholded_tmap, adjacency_matrix)

    [n_verts, n_timepoints] = size(thresholded_tmap);

	% To begin with we don't look for temporal contiguity.  We label
    % spatially contiguous clusters in each timepoint with a unique
    % label.
    spatial_cluster_labels = zeros(n_verts.(chi), n_timepoints);
    
    running_cluster_count = 0;
    for t = 1:n_timepoints
        
        % We're interested in identifying contiguous clusters, so we'll
        % forget adjacency information for sub-threshold vertices.
        masked_adjacency_matrix = adjacency_matrix;
        masked_adjacency_matrix(thresholded_tmap(:, t) == 0, :                         ) = 0;
        masked_adjacency_matrix(                         :, thresholded_tmap(:, t) == 0) = 0;
        
        % A cell array of components at this timepoint. Each component is a vector of
        % vertex names.
        component_list_this_t = connected_components(masked_adjacency_matrix);
        n_clusters_this_t = numel(component_list_this_t);
        
        for component_i = 1:n_clusters_this_t
            spatial_cluster_labels(component_list_this_t(component_i), t) = running_cluster_count + component_i;
        end
        
        running_cluster_count = running_cluster_count + n_clusters_this_t;
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
            for overlap_vi = vis_overlap
               cluster_id_pairs_to_merge = [ ...
                   cluster_id_pairs_to_merge; ...
                   [spatial_cluster_labels(overlap_vi, t), spatial_cluster_labels(overlap_vi, t+1)]];
            end
            cluster_id_pairs_to_merge = unique(cluster_id_pairs_to_merge);
            
            if numel(cluster_id_pairs_to_merge) > 0
                merging_done_last_pass = true;
            end

            % merge clusters
            for cluster_pair_i = 1:size(cluster_id_pairs_to_merge, 2);
                cluster_id_pair = cluster_id_pairs_to_merge(:, cluster_pair_i);

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
            for overlap_vi = vis_overlap
               cluster_id_pairs_to_merge = [ ...
                   cluster_id_pairs_to_merge; ...
                   [spatial_cluster_labels(overlap_vi, t), spatial_cluster_labels(overlap_vi, t-1)]];
            end
            cluster_id_pairs_to_merge = unique(cluster_id_pairs_to_merge);
            
            if numel(cluster_id_pairs_to_merge) > 0
                merging_done_last_pass = true;
            end

            % merge clusters
            for cluster_pair_i = 1:size(cluster_id_pairs_to_merge, 2);
                cluster_id_pair = cluster_id_pairs_to_merge(:, cluster_pair_i);

                % relabel all vertices of the matching clusters in both
                % layers to have the same label
                spatial_cluster_labels(spatial_cluster_labels == max(cluster_id_pair)) = min(cluster_id_pair);

                % relabel all to-merge clusters too
                cluster_id_pairs_to_merge(cluster_id_pairs_to_merge == max(cluster_id_pair)) = min(cluster_id_pair);
            end
        end%for
    end%while
    
    % Relabel clusters
    
    remaining_cluster_labels = unique(spatial_cluster_labels);
    % these will be sorted
    
    n_clusters = numel(remaining_cluster_labels);
    
    % renumber clusters in ascending order to prevent collisions
    for rem_cluster_lab_i = 1:n_clusters
        rem_cluster_lab = remaining_cluster_labels(rem_cluster_lab_i);
        spatial_cluster_labels(spatial_cluster_labels == rem_cluster_lab) = rem_cluster_lab_i;
    end
end

