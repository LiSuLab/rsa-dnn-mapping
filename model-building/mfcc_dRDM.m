function dRDM = mfcc_dRDM(distance_type, frame_cap);

    %% Paths

    % Change these values
    input_dir = '/Users/cai/Desktop/scratch/py_out/cepstral-coefficients';


    %% Load in the features

    
    N_coeffs_per_type = 12;
    for coeff_type = 'CDA'
        for coeff_i_within_type = 1:N_coeffs_per_type
            
            coeff_name = sprintf('%s%02d', coeff_type, coeff_i_within_type);
            
            file_name = sprintf('cepstral-coefficients-%s.mat', coeff_name);
            file_path = fullfile(input_dir, file_name);
            
            coeffs.(coeff_name) = load(file_path);
            
        end 
    end
    
    word_list = fieldnames(coeffs.C01);
    % CONDITIONS ARE IN ALPHABETICAL ORDER OF WORDS
    word_list = sort(word_list);
    n_words = numel(word_list);    
    n_frames = numel(coeffs.C01.(word_list{1}));


    %% Produce RDMs
    
    coeff_types_for_model = 'C';
    coeff_indices_for_model = 1:12;
    n_coeffs_for_model = numel(coeff_types_for_model) * numel(coeff_indices_for_model);

    % Frames per window. For HTK, each frame is 10ms.
    window_width_in_frames = 2;

    n_windows = n_frames - window_width_in_frames + 1;

    % We scan over the coefficients with a sliding window
    for frame_i = 1:min([n_windows, frame_cap])
        
        % The window of indices relative to the coefficient timelines
        window = frame_i:frame_i+window_width_in_frames-1;
        
        % preallocate data for this frame
        % should eventually be a words-by-(window_width*n_coeffs) data matrix
        data_this_frame = nan(n_words, window_width_in_frames * n_coeffs_for_model);
        for word_i = 1:n_words
            word = word_list{word_i};
            
            data_this_word = [];
            
            overall_coeff_i = 1;
            for coeff_type = coeff_types_for_model
                for coeff_i_within_type = coeff_indices_for_model
                    coeff_name = sprintf('%s%02d', coeff_type, coeff_i_within_type);
                    data_this_word = [ ...
                        data_this_word, ...
                        coeffs.(coeff_name).(word)(window)];
                    overall_coeff_i = overall_coeff_i + 1;
                end
            end
            
            data_this_frame(word_i, :) = data_this_word;
        end
        
        RDM_this_frame = pdist(data_this_frame, distance_type);
        
        dRDM(frame_i).RDM  = RDM_this_frame;
        dRDM(frame_i).Name = sprintf('mfcc-%02d', frame_i);

    end%for

end
