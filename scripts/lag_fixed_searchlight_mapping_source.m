function [mapsPath] = lag_fixed_searchlight_mapping_source(subjectName, chi, RDMPath, slMask, modelRDM, model_lag_ms, STCMetadatas, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.util.*

    %% File paths

    mapsDir = fullfile(userOptions.rootPath, 'Maps');
    mapsFileName = fprintf('%s_rMesh_%s_lagfix%dms_sub%s-%sh.stc', ...
        userOptions.analysisName, ...
        spacesToUnderscores(modelRDM.name), ...
        model_lag_ms, ...
        subjectName, ...
        lower(chi));
    mapsPath = fullfile(mapsDir, mapsFileName);

    promptOptions.functionCaller = 'searchlightMapping_source';
    promptOptions.defaultResponse = 'S';
    promptOptions.checkFiles(1).address = mapsPath;

    overwriteFlag = overwritePrompt(userOptions, promptOptions);

    %% Apply searchlight

    if overwriteFlag

        [slSpecs, slSTCMetadatas] = getSearchlightSpec(STCMetadatas, userOptions);

        searchlightRDMs = directLoad(RDMPath, 'searchlightRDMs');

        modelRDM_utv = squeeze(unwrapRDMs(vectorizeRDMs(modelRDM)));
        
        n_mask_vertices = numel(indexMask.vertices);

        % The number of positions the sliding window will take.
        nWindowPositions = size(slSpecs.windowPositions, 1);

        %% map the volume with the searchlight

        % Preallocate looped matrices for speed
        thisSubjectRs = zeros(n_mask_vertices, nWindowPositions);

        % For display purposes
        nVertsSearched = 0;

        % Search the vertices
        for v_i = 1:n_mask_vertices

            % Search through time
            window_i = 0;
            for window = slSpecs.windowPositions'
                
                % This counts the windows (1-indexed) as they are
                % considered.
                window_i = window_i + 1;

                % Get this RDM by vertex and window indices, as that's how 
                % it was stored.
                patchRDM = searchlightRDMs(v_i, window_i).RDM;

                rs = corr(patchRDM', modelRDM_utv', 'type', userOptions.RDMCorrelationType, 'rows', 'pairwise');

                % Store results to be retured.
                thisSubjectRs(v_i, window_i) = rs;

            end%for:window

            % Indicate progress every once in a while...
            nVertsSearched = nVertsSearched + 1;
            if mod(nVertsSearched, 200) == 0
                prints('%d vertices searched: %d%% complete', nVertsSearched, floor((nVertsSearched / numel(indexMask.vertices)) * 100));
            end%if

        end%for:v

        % Wrap in a metadata struct
        
        rSTCStruct          = slSTCMetadatas.(chi);
        rSTCStruct.vertices = slMask.vertices;
        % thisSubjectRs contains only the data inside the mask, but since the
        % vertices are stored in this struct, that should be ok.
        rSTCStruct.data     = thisSubjectRs(:,:);
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
