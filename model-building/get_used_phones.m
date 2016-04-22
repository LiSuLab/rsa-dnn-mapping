% Gets the list of all used phones from the list of segmentations.
% Returns a cell array phones.
function phones = get_used_phones(segmentations)
    
    words = fieldnames(segmentations);
    n_words = numel(words);
    
    % list of phones
    phones = {};
    
    % Get all phones used with duplicates
    for word_i = 1:n_words
        word = words{word_i};
        
        n_segments = numel(segmentations.(word));
        for segment_i = 1:n_segments
            segment = segmentations.(word)(segment_i);
            
            phones = [phones {segment.label}];
        end
    end
    
    % Remove duplicates
    phones = unique(phones);
    
    % Remove 'sil' and 'sp', if they crept in
    permitted_phone_is = cellfun( ...
        @(x)( ...
            ~strcmpi(x, 'sil') ...
            && ~strcmpi(x, 'sp') ...
        ), ...
        phones);
    phones = phones(permitted_phone_is);
    
    % Sort
    phones = sort(phones);
        
end
