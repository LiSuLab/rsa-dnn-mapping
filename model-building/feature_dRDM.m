function dRDM = feature_dRDM(distance_type, frame_cap)

    if ~exist('distance_type', 'var'), distance_type = 'Correlation'; end
    
    feature_traces = feature_trace_extraction();
    
    features = fieldnames(feature_traces);
    n_features = numel(features);
    
    [n_words, n_timepoints] = size(feature_traces.(features{1}));
    
    dRDM = struct();
    for t = 1:min([n_timepoints, frame_cap])
    
        % For each timepoint, build a word-by-nfeatures data matrix
        data_this_t = nan(n_words, n_features);
        for feature_i = 1:n_features
            feature = features{feature_i};
            for word_i = 1:n_words
                data_this_t(word_i, feature_i) = feature_traces.(feature)(word_i, t);
            end
        end
        
        RDM_this_t = pdist(data_this_t, distance_type);
        dRDM(t).RDM  = RDM_this_t;
        dRDM(t).Name = sprintf('feature_%02d', frame_i);
        
    end
end