function [dRDMs, cm] = compare_dRDMs();

    rdm_corr_type = 'Spearman';
    rdm_type      = 'correlation';
    frame_cap = 25;
    
    % Build model struct
    rdm_i = 1;
    dRDMs(:, rdm_i) = mfcc_dRDM(                         rdm_type, frame_cap); rdm_i = rdm_i + 1;
    dRDMs(:, rdm_i) = dynamic_hidden_layer_models('2',   rdm_type, frame_cap); rdm_i = rdm_i + 1;
    dRDMs(:, rdm_i) = dynamic_hidden_layer_models('3',   rdm_type, frame_cap); rdm_i = rdm_i + 1;
    dRDMs(:, rdm_i) = dynamic_hidden_layer_models('4',   rdm_type, frame_cap); rdm_i = rdm_i + 1;
    dRDMs(:, rdm_i) = dynamic_hidden_layer_models('5',   rdm_type, frame_cap); rdm_i = rdm_i + 1;
    dRDMs(:, rdm_i) = dynamic_hidden_layer_models('6',   rdm_type, frame_cap); rdm_i = rdm_i + 1;
    dRDMs(:, rdm_i) = dynamic_hidden_layer_models('7BN', rdm_type, frame_cap); rdm_i = rdm_i + 1;
    dRDMs(:, rdm_i) = triphone_dRDM(                     rdm_type, frame_cap); rdm_i = rdm_i + 1;
    dRDMs(:, rdm_i) = phone_dRDM(                        rdm_type, frame_cap); rdm_i = rdm_i + 1;
    dRDMs(:, rdm_i) = feature_dRDM(                      rdm_type, frame_cap); rdm_i = rdm_i + 1;
    
    [n_frames, n_models] = size(dRDMs);
    model_size = numel(dRDMs(1).RDM);
    
    for t = 1:n_frames
        
        % Get RDMs in for this frame
        models_this_frame = nan(n_models, model_size);
        for m = 1:n_models
            models_this_frame(m, :) = dRDMs(t, m).RDM;
        end
        
        % dynamic correlation matrix
        [cm(:, :, t), pm(:, :)] = corr( ...
            models_this_frame', ...
            'type', rdm_corr_type);
        
        % dynamic distance matrix
        dm = 1 - cm(:, :, t);
        % where there are nans, put large distance
        dm(isnan(dm)) = 1;
        % make sure every point is at distance zero from itself
        dm(find(eye(size(dm, 1)))) = 0; %#ok<FNDSB> % 'fixing' this actually breaks it
        
        figure;
        subplot(3,1,1);
        
        points(:, :) = mdscale( ...
             dm(:, :), 2 ...
            ,'Criterion', 'sammon');
        scatter(points(:, 1), points(:, 2), 140, 'filled');
        % Label them
        for m = 1:n_models
            text(points(m, 1), points(m, 2), dRDMs(t, m).Name);
        end
        
        subplot(3, 1, 2);
        imagesc(cm(:,:, t)); colorbar;
        
        subplot(3, 1, 3);
        imagesc(pm < 0.05);
        
    end
    
end
