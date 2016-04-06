% Given an adjacency matrix, some spatiotemporal maps, and a
% cluster-forming threshold, this will return labelled spatiotemporal
% clusters.
function [labelled_spatiotemporal_clusters, vertex_level_threshold] = identify_spatiotemporal_clusters(adjacency_matrix, group_maps, cluster_forming_threshold)
    vertex_level_threshold = quantile(group_maps(:), 1-cluster_forming_threshold);
    thresholded_map = (group_maps > vertex_level_threshold);
    labelled_spatiotemporal_clusters = label_spatiotemporal_clusters(thresholded_map, adjacency_matrix);
end

%% SUBFUNCTIONS

% Given a binary spatiotemporal map and an adjacency matrix, this will
% return a map which looks the same as the binary map, but where each
% spatiotemporally contiguous cluster has a unique integer label.
function spatial_cluster_labels = label_spatiotemporal_clusters(binary_map, adjacency_matrix)

    [n_verts, n_timepoints] = size(binary_map);
    
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
        super_threshold_vs_iwm = find(binary_map(:, t));
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
