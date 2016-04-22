function dRDM = phone_dRDM(distance_type, frame_cap)
    
    phone_models = phone_activations();
    phone_models = orderfields(phone_models);
    
    phones = fieldnames(phone_models);
    n_phones = numel(phones);
    
    [n_words, n_timepoints] = size(phone_models.(phones{1}));
    
    dRDM = struct();
    for t = 1:min([n_timepoints, frame_cap])
    
        % For each timepoint, build a word-by-nfeatures data matrix
        data_this_t = nan(n_words, n_phones);
        for phone_i = 1:n_phones
            phone = phones{phone_i};
            for word_i = 1:n_words
                data_this_t(word_i, phone_i) = phone_models.(phone)(word_i, t);
            end
        end
        
        RDM_this_t = pdist(data_this_t, distance_type);
        dRDM(t).RDM  = RDM_this_t;
        dRDM(t).Name = sprintf('phone-%02d', t);
        
    end
end
