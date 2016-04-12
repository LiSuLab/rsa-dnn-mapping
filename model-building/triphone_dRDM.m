function dRDM = triphone_dRDM(distance_type)

    if ~exist('distance_type', 'var'), distance_type = 'Correlation'; end
    
    triphone_data = load('/Users/cai/Desktop/scratch/py_out/triphone_likelihoods/actual_triphone_values.mat');
    triphone_data = orderfields(triphone_data);
    
    words = fieldnames(triphone_data);
    n_words = numel(words);
    
    % This data starts from frame 2, so we must add an additional blank
    % frame to begin with.  This is just an aspect of using HTK's output.
    [n_timepoints_trimmed, n_triphones] = size(triphone_data.(words{1}));
    
    dRDM = struct();
    
    dRDM(1).RDM = squareform(zeros(n_words, n_words));
    
    for t = 1:n_timepoints_trimmed
       
        data_this_frame = nan(n_words, n_triphones);
        for word_i = 1:n_words
            word = words{word_i};
            data_this_frame(word_i, :) = triphone_data.(word)(t, :);
        end
        
        dRDM(t+1).RDM = pdist(data_this_frame, distance_type);
        
    end

end
