% Given an adjacency matrix, some spatiotemporal maps, and a
% cluster-forming threshold, this will return labelled spatiotemporal
% clusters.
function [labelled_spatiotemporal_clusters, vertex_level_threshold] = identify_spatiotemporal_clusters(adjacency_matrix, group_maps, cluster_forming_threshold)
    vertex_level_threshold = quantile(group_maps(:), 1-cluster_forming_threshold);
    thresholded_map = (group_maps > vertex_level_threshold);
    labelled_spatiotemporal_clusters = label_spatiotemporal_clusters(thresholded_map, adjacency_matrix);
end

