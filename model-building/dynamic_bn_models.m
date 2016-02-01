function dRDM = dynamic_bn_models()

    bn_activations = load('/Users/cai/Desktop/scratch/py_out/hidden_layer_7BN_activations');
    
    words = fieldnames(bn_activations);
    n_words = numel(words);
    
    shortest_word_length = inf;
    for word_i = 1:n_words
        word = words{word_i};
        [word_length, n_bn_nodes] = size(bn_activations.(word));
        shortest_word_length = min(shortest_word_length, word_length);
    end
    
    dRDM = struct();
    for t = 1:shortest_word_length
       data_this_timepoint = nan(n_words, n_bn_nodes);
       for word_i = 1:n_words
           word = words{word_i};
           word_activation = bn_activations.(word);
           data_this_timepoint(word_i, :) = word_activation(t, :);
       end
       RDM_this_timepoint = pdist(data_this_timepoint, 'correlation');
       dRDM(t).RDM = RDM_this_timepoint;
    end

end%function
