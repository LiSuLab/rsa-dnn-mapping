function [dCM, dP, dPoints] = compare_dRDMs();

    rdm_corr_type = 'Spearman';
    rdm_type      = 'correlation';

    phone_RDMs    = phone_dRDM(rdm_type);
    feature_RDMs  = feature_dRDM(rdm_type);
    triphone_RDMs = triphone_dRDM(rdm_type);
    
    % Build model struct
    
    % Start with hidden layers
    dRDMs = dnn_layer_dRDMs(rdm_type);
    
    n_frames = numel(dRDMs);
    
    for t = 1:n_frames
        dRDMs(t).Phones.Name = sprintf('phones,t=%02d', t);
        dRDMs(t).Phones.RDM = phone_RDMs(t).RDM;
        
        dRDMs(t).Features.Name = sprintf('features,t=%02d', t);
        dRDMs(t).Features.RDM = feature_RDMs(t).RDM;
        
        dRDMs(t).Triphones.Name = sprintf('triphones,t=%02d', t);
        dRDMs(t).Triphones.RDM = triphone_RDMs(t).RDM;
    end
    
    model_names = fieldnames(dRDMs);
    n_models = numel(model_names);
    
    model_size = numel(dRDMs(1).(model_names{1}).RDM);
    
    c = [ ...
        0., 0., 0.; ...
        .2, 0., 0.; ...
        .4, 0., 0.; ...
        .6, 0., 0.; ...
        .8, 0., 0.; ...
        1., 0., 0.; ...
        0., 1., 0.; ...
        0., 1., 1.; ...
        0., 0., 1.];
    
    all_models = [];
    for t = 1:n_frames
        
        % Get RDMs in for this frame
        models_this_frame = nan(n_models, model_size);
        for m = 1:n_models
            model_name = model_names{m};
            models_this_frame(m, :) = dRDMs(t).(model_name).RDM;
            all_models = [all_models; models_this_frame(m, :)];
        end
        
        % dynamic correlation matrix
        [dCM(:, :, t), dP(:, :, t)] = corr( ...
            models_this_frame', ...
            'type', rdm_corr_type);
        
        % dynamic distance matrix
        dm = 1 - dCM(:, :, t);
        % where there are nans, put large distance
        dm(isnan(dm)) = 1;
        % make sure every point is at distance zero from itself
        dm(find(eye(size(dm, 1)))) = 0; %#ok<FNDSB> % 'fixing' this actually breaks it
        dDM(:, :, t) = dm;
        
        % Look at this gross abuse of the lack of static typing in Matlab
        if t == 1
            start = 'random';
        else
            start = dPoints(:, :, t-1);
        end
        
        dPoints(:, :, t) = mdscale( ...
            dDM(:, :, t), 2 ...
            ,'Criterion', 'sammon' ...
            ...%,'Weights', 0.5*blkdiag(ones(6),ones(3)) + 0.5*ones(9) ...
            ...%,'Options', struct( ...
            ...%    'MaxIter', 300) ...
            ...%,'Start', start ...
            );
        
        figure;
        subplot(2,1,1);
        scatter(dPoints(:, 1, t), dPoints(:, 2, t), 140, c, 'filled');
        % Label them
        for m = 1:n_models
            text(dPoints(m,1,t), dPoints(m,2,t), dRDMs(t).(model_names{m}).Name);
        end
        subplot(2,1,2);
        imagesc(dCM(:,:,t)); colorbar;
        
    end
    
    % All Frames Together
    [CM, P] = corr(all_models', 'type', rdm_corr_type);
    dm = 1-CM;
    dm(isnan(dm)) = 1;
    dm(find(eye(size(dm, 1)))) = 0; %#ok<FNDSB>
    figure;
    subplot(2,1,1);
    m = mdscale(dm, 2, 'Criterion', 'sammon', 'Start', 'random');
    scatter(m(:, 1), m(:, 2), 140, repmat(c, n_frames, 1), 'filled');
    subplot(2,1,2);
    imagesc(dm); colorbar;
    
end
