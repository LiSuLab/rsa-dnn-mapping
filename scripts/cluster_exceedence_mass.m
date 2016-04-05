% Given some spatiotemporal maps and some identically-sized cluster labels,
% this will return a cluster-label-indexed vector of cluster exceedence
% masses.
function cluster_stats = cluster_exceedence_mass(labelled_spatiotemporal_clusters, group_maps)
    
    cluster_ids = unique(labelled_spatiotemporal_clusters);
    % The cluster whose id is zero is not a cluster at all, so we delete it
    % here. It's the background!
    cluster_ids = cluster_ids(cluster_ids > 0);
    
    % preallocate
    cluster_stats = nan(size(cluster_ids));
    
    for cluster_i = cluster_ids'
        % cluster exceedence mass
        this_cluster_location = (labelled_spatiotemporal_clusters == cluster_i);
        cluster_exceedences = group_maps(this_cluster_location) - vertex_level_threshold;
        cluster_stats(cluster_i) = sum(cluster_exceedences(:));
    end
end
