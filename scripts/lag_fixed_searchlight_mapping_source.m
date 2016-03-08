function [mapsPath] = lag_fixed_searchlight_mapping_source(chi, data_RDM_paths, name_prefix, slMask, modelRDM, model_lag_ms, STCMetadatas, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.util.*

    %% File paths

    mapsDir = fullfile(userOptions.rootPath, 'Maps');
    mapsFileName = sprintf('%s_lagfix_%dms-%sh.stc', ...
        name_prefix, ...
        model_lag_ms, ...
        lower(chi));
    mapsPath = fullfile(mapsDir, mapsFileName);

    promptOptions.functionCaller = 'lag_fixed_searchlight_mapping_source';
    promptOptions.defaultResponse = 'S';
    promptOptions.checkFiles(1).address = mapsPath;

    overwriteFlag = overwritePrompt(userOptions, promptOptions);

    %% Apply searchlight

    if overwriteFlag

        [slSpecs, slSTCMetadatas] = getSearchlightSpec(STCMetadatas, userOptions);

        modelRDM_utv = squeeze(unwrapRDMs(vectorizeRDMs(modelRDM)));
        
        n_mask_vertices = numel(slMask.vertices);

        % The number of positions the sliding window will take.
        nWindowPositions = size(slSpecs.(chi).windowPositions, 1);

        %% map the volume with the searchlight

        % Preallocate looped matrices for speed
        % Preallocate as zeros as it will be saved as an STC.
        r_mesh = zeros(n_mask_vertices, nWindowPositions);

        % Search through time
        parfor window_i = 1:nWindowPositions
            prints('\tWindow %d/%d', window_i, nWindowPositions);

            searchlightRDMs = directLoad(data_RDM_paths(window_i).(chi));

            % Search the vertices
            for v_i = 1:n_mask_vertices

                % Store results to be retured.
                r_mesh(v_i, window_i) = corr( ...
                    ...% Get this RDM by vertex and window indices, as 
                    ...% that's how it was stored.
                    searchlightRDMs(v_i, :)', ...
                    modelRDM_utv', ...
                    'type', userOptions.RDMCorrelationType, ...
                    'rows', 'pairwise');

            end%for:v_i
        end%for:window

        
        %% Wrap in a metadata struct
        
        rSTCStruct          = slSTCMetadatas.(chi);
        rSTCStruct.vertices = slMask.vertices;
        % r_mesh contains only the data inside the mask, but since the
        % vertices are stored in this struct, that should be ok.
        rSTCStruct.data     = r_mesh;
        % Correct lag.
        % In the case where model_lag is greater than zero, the model is
        % making predictions about the data at a positively shifted lag.
        % Therefore we must negatively shift the lag in the result - e.g.
        % if we have a model for +10ms audio, the matching this with +0ms
        % data is like making a prediction before the audio was heard:
        % -10ms.
        rSTCStruct.tmin     = rSTCStruct.tmin - (model_lag_ms / 1000);
        rSTCStruct.tmax     = rSTCStruct.tmax - (model_lag_ms / 1000);
        

        %% Saving r-maps and RDM maps

        prints('Writing r-map %s.', mapsPath);
        mne_write_stc_file1(mapsPath, rSTCStruct);
        
    else
        prints('Searchlight already applied, skipping it.');
    end

end%function
