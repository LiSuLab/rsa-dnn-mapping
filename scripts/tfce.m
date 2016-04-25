% Applies threshold-free cluster enhancement to stc data (method of Smith
% and Nichols Neuroimage 2009).
%
% Original author: Cai Wingfield 2016-04
function map_tfce = tfce(map_raw, adjacency_matrix_iwm)

    import rsa.*
    import rsa.meg.*
    import rsa.util.*
    
    E = 0.5;
    H = 2;
    
    % TODO:
    
    % change in height
    dh = 0.1;
    
    % The maximum height
    max_h = max(map_raw(:));
    
    heights = 0 : dh : max_h;
    n_heights = numel(heights);
    
    [n_vertices, n_timepoints] = size(map_raw);
    
    % Precompute cluster extents at each height, at each datapoint
    
    % this will hold, for each vertex p, the tfce value
    %   \int_{h=0}^{h=h_p} extent_{cluster_p}(h)^E h^H dh
    % (Smith & Nichols 2009)
    % as approximated by the descretised sum
    %   \sum_k (extent_{cluster_p}(k dh))^E (k dh)^H dh
    % since h = k*dh for appropriately small finite dh, index k
    
    map_tfce = zeros(n_vertices, n_timepoints);
    
    for h_i = 1:n_heights
        
        % this h threshokld
        h = heights(h_i);
      
        % the map thresholded a this h
        h_thresholded_map = map_raw > h;
        
        % label the clusters
        h_labelled_clusters = label_spatiotemporal_clusters( ...
            h_thresholded_map, ...
            adjacency_matrix_iwm);
        
        % list of cluster extents, by cluster label
        h_cluster_extents = cluster_extent(h_labelled_clusters);
        
        % a map of clusters with in-cluster vertices labelled with their
        % extents
        h_extent_labelled_clusters = zeros(n_vertices, n_timepoints);
        for cluster_i = 1:numel(h_cluster_extents)
            h_extent_labelled_clusters(h_labelled_clusters == cluster_i) = h_cluster_extents(cluster_i);
        end
        
        summand_map = (h_extent_labelled_clusters .^ E) .* (h^H) .* dh;
        
        map_tfce = map_tfce + summand_map;
        
    end

end%function