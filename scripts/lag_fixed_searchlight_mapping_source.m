function [mapsPath] = lag_fixed_searchlight_mapping_source(chi, data_RDM_paths, slMask, modelRDM, model_lag_ms, STCMetadatas, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.util.*

    %% File paths

    mapsDir = fullfile(userOptions.rootPath, 'Maps');
    mapsFileName = sprintf('%s_rMesh_lagfix_%dms-%sh.stc', ...
        userOptions.analysisName, ...
        model_lag_ms, ...
        lower(chi));
    mapsPath = fullfile(mapsDir, mapsFileName);

    promptOptions.functionCaller = 'searchlightMapping_source';
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
        % Preallocate as zeros asthis will be saved as an STC.
        r_mesh = zeros(n_mask_vertices, nWindowPositions);

        % Search through time
        window_i = 0;
        for window = slSpecs.(chi).windowPositions'
                
            % This counts the windows (1-indexed) as they are
            % considered.
            window_i = window_i + 1;
            
            prints('Window %d', window_i);

            searchlightRDMs = directLoad(data_RDM_paths(window_i).(chi));

            % For display purposes
            nVertsSearched = 0;

            % Search the vertices
            for v_i = 1:n_mask_vertices

                % Get this RDM by vertex and window indices, as that's how 
                % it was stored.
                patchRDM = squeeze(searchlightRDMs(v_i, :));

                rs = corr(patchRDM', modelRDM_utv', 'type', userOptions.RDMCorrelationType, 'rows', 'pairwise');

                % Store results to be retured.
                r_mesh(v_i, window_i) = rs;

                % Indicate progress every once in a while...
                nVertsSearched = nVertsSearched + 1;
                if mod(nVertsSearched, 200) == 0
                    prints('%d vertices searched, %d percent complete', nVertsSearched, floor((nVertsSearched / numel(slMask.vertices)) * 100));
                end%if

            end%for:v_i

        end%for:window

        % Wrap in a metadata struct
        
        rSTCStruct          = slSTCMetadatas.(chi);
        rSTCStruct.vertices = slMask.vertices;
        % r_mesh contains only the data inside the mask, but since the
        % vertices are stored in this struct, that should be ok.
        rSTCStruct.data     = r_mesh(:,:);
        % Correct for lag.
        rSTCStruct.tmin     = rSTCStruct.tmin - model_lag_ms;

        %% Saving r-maps and RDM maps

        prints('Writing r-map %s.', mapsPath);
        gotoDir(mapsDir);
        mne_write_stc_file1(mapsPath, rSTCStruct);
        
    else
        prints('Searchlight already applied, skipping it.');
    end

end%function
