function dRDM = dynamic_hidden_layer_models(layer_name, distance_type, frame_cap)

    bn_activations = load(sprintf('/Users/cai/Desktop/scratch/py_out/hidden_layer_%s_activations.mat', layer_name));
    
    words = fieldnames(bn_activations);
    % CONDITIONS ARE IN ALPHABETICAL ORDER OF WORDS
    words = sort(words);
    n_words = numel(words);
    
    shortest_word_length = inf;
    for word_i = 1:n_words
        word = words{word_i};
        [word_length, n_bn_nodes] = size(bn_activations.(word));
        shortest_word_length = min(shortest_word_length, word_length);
    end
    
    dRDM = struct();
    for t = 1:min([shortest_word_length, frame_cap])
       data_this_timepoint = nan(n_words, n_bn_nodes);
       for word_i = 1:n_words
           word = words{word_i};
           word_activation = bn_activations.(word);
           data_this_timepoint(word_i, :) = word_activation(t, :);
       end
       RDM_this_timepoint = pdist(data_this_timepoint, distance_type);
       dRDM(t).RDM  = RDM_this_timepoint;
       dRDM(t).Name = sprintf('%s-%02d', layer_name, t);
    end

end%function
