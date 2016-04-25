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
