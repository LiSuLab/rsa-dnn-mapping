% Produces sliding-window-matched plots for the instances of each phone for
% each word.
function [phone_models] = phone_activations()
    
    load_dir = '/imaging/cw04/CSLB/Lexpro/Analysis_DNN/Models/';

    segmentations = load(fullfile(load_dir, 'triphone_boundaries.mat'));
    segmentations = orderfields(segmentations);

    bn26 = load(fullfile(load_dir, 'hidden_layer_7BN_activations.mat'));
    bn26 = orderfields(bn26);
    
    phones = get_used_phones(segmentations);
    words = fieldnames(segmentations);
    % CONDITIONS ARE IN ALPHABETICAL ORDER OF WORDS
    words = sort(words);
    
    n_phones = numel(phones);
    n_words = numel(words);
    
    % These phone plots will model the node plots, so we will produce
    % word-by-time plots for each phone.
    
    max_word_length = 1;
    for word_i = 1:n_words
       word = words{word_i};
       max_word_length = max(max_word_length, size(bn26.(word), 1));
    end
    
    phone_models = struct();
    
    for phone_i = 1:n_phones
        phone = phones{phone_i};
        
        % There are some phones we don't care about for now
        if strcmpi(phone, 'sil') || strcmpi(phone, 'sp')
            continue;
        end
        
        % Initialise empty matrix
        phone_models.(phone) = nan(n_words, max_word_length);
        
        for word_i = 1:n_words
            word = words{word_i};
            this_word_length = size(bn26.(word), 1);
            
            % move a virtual sliding window throughout the segmentatino of
            % the word to get a model activation map for this phone.
            
            % Most words don't have most phones, so we can quickly check
            % that and just insert all zeros where necessary.
            
            this_word_phones = unique({ segmentations.(word)(:).label });
            
            % Initialise with zeros to begin with
            phone_profile = zeros(this_word_length, 1);
            
            if any(strcmpi(this_word_phones, phone))
                % The phone is in the word
                
                % Get segmentations indices whose label is this phone
                phone_segment_is = find(strcmpi(phone, {segmentations.(word)(:).label}));
                
                % For each time this phone exists in this word
                for phone_segment_i = phone_segment_is
                    phone_segment_100ns = [ ...
                        segmentations.(word)(phone_segment_i).onset, ...
                        segmentations.(word)(phone_segment_i).offset ...
                    ];
                    % Convert out of int64
                    phone_segment_100ns = double(phone_segment_100ns);
                    
                    for frame_i = 1:this_word_length
                        frame_segment = frame2100ns(frame_i);
                        % Determine to what degree this frame overlaps with
                        % the phone
                        overlap_this_frame = get_interval_overlap(frame_segment, phone_segment_100ns);
                        % Write this into the profile.
                        phone_profile(frame_i) = phone_profile(frame_i) + overlap_this_frame;
                    end
                end
            end
            
            % Insert the profile into the big struct
            phone_models.(phone)(word_i, 1:this_word_length) = phone_profile';
            
        end
    end
end%function

%% SUB FUNCTIONS

% converts miliseconds into the units used by HTK to perform segmentations
% (which are each 100ns)
function v_100ns = ms2100ns(v_ms)
    MS_2_100NS = 10000; % how many 100ns are there in a ms
    v_100ns = v_ms * MS_2_100NS;
end%function

% converts frames into the units used by HTK to perform segmentations
% (which are each 100ns)
function seg_100ns = frame2100ns(f)
    FRAME_STEP_MS  = 10;
    FRAME_WIDTH_MS = 25;
    % convert to miliseconds (the -1 is from the 1-indexing of frames vs
    % 0-indexing of ms)
    seg_ms = (f - 1) * FRAME_STEP_MS;
    seg_ms = [seg_ms, seg_ms + FRAME_WIDTH_MS];
    seg_100ns = ms2100ns(seg_ms);
end%function

