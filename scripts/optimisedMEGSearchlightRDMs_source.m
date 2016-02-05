function [RDMsPaths, slSTCMetadatas, slSpecs] = optimisedMEGSearchlightRDMs_source(meshPaths, slMasks, adjacencyMatrix, STCMetadatas, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.util.*

    %% Constants

    usingMasks = ~isempty(userOptions.maskNames);

    nSubjects = numel(userOptions.subjectNames);

    %% File paths

    RDMsDir = fullfile(userOptions.rootPath, 'RDMs');

    file_i = 1;
    for subject_i = 1:nSubjects
        thisSubjectName = userOptions.subjectNames{subject_i};
        for chi = 'LR'
            if usingMasks
                RDMsFile = ['searchlightRDMs_masked_', thisSubjectName, '-' lower(chi) 'h.mat'];
            else
                RDMsFile = ['searchlightRDMs_',        thisSubjectName, '-' lower(chi) 'h.mat'];
            end
            RDMsPaths(subject_i).(chi) = fullfile(RDMsDir, RDMsFile);

            % We'll check all the files to be saved to see if they have already
            % been saved.
            promptOptions.checkFiles(file_i).address = RDMsPaths(subject_i).(chi);
            file_i = file_i + 1;
        end
    end

    promptOptions.functionCaller = 'MEGSearchlightRDMs_source';
    promptOptions.defaultResponse = 'S';

    overwriteFlag = overwritePrompt(userOptions, promptOptions);


    %% Apply searchlight

    [slSpecs, slSTCMetadatas] = getSearchlightSpec(STCMetadatas, userOptions);

    parfor subject_i = 1:nSubjects
        thisSubjectName = userOptions.subjectNames{subject_i};

        % Work on each hemisphere separately
        for chi = 'LR'

            % We'll only do the searchlight if we haven't already done so,
            % unless we're being told to overwrite.
            if exist(RDMsPaths(subject_i).(chi), 'file') && ~overwriteFlag
                prints('Searchlight already performed in %sh hemisphere of subject %d. Skipping.', lower(chi), subject_i);
                
            else
                prints('Shining RSA searchlight in the %sh source mesh of subject %d of %d (%s)...', lower(chi), subject_i, nSubjects, thisSubjectName);

                single_hemisphere_searchlight( ...
                    slSpecs.(chi), ...
                    slMasks([slMasks.chi] == chi), ...
                    meshPaths(subject_i).(chi), ...
                    RDMsPaths(subject_i).(chi), ...
                    adjacencyMatrix, ...
                    userOptions);

                %% Done
                prints('Done with subeject %d''s %sh side.', subject_i, lower(chi));

            end%if:overwrite
        end%for:chi
    end%for:subject

end%function

%%%%%%%%%%%%%%%%%%
%% Subfunctions %%
%%%%%%%%%%%%%%%%%%

% Computes and saves searchlight RDMs for a single hemisphere of a single
% subject.
function single_hemisphere_searchlight(slSpec, indexMask, meshPath, RDMsPath, adjacencyMatrix, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.util.*

    maskedMeshes = directLoad(meshPath, 'sourceMeshes');

    [nVertices_masked, nTimePoints_data, nConditions, nSessions] = size(maskedMeshes);
    
    n_rdm_entries = numel(vectorizeRDM(zeros(nConditions)));

    % The number of positions the sliding window will take.
    nWindowPositions = size(slSpec.windowPositions, 1);

    %% map the volume with the searchlight

    % Preallocate looped matrices for speed
    searchlightRDMs = nan(nVertices_masked, nWindowPositions, n_rdm_entries);

    % For display purposes only
    nVertsSearched = 0;

    % Search the vertices
    for v_i = 1:numel(indexMask.vertices)
        % v_i loops through the *indices of* vertices in the mask
        % v is the vertex itself

        v = indexMask.vertices(v_i);

        % Determine which vertexes are within the radius of the currently-picked vertex
        searchlightPatchVs = [v, adjacencyMatrix(v,:)];

        % Restrict to verticies inside mask.
        % This also removes any nans.
        searchlightPatchVs = intersect(searchlightPatchVs, indexMask.vertices);

        % Now we need to convert the vertices into vertex *indices*.  For
        % example, our mask may be vertices [5, 6, 7], but since there will
        % only be three datapoints inside each of the masked meshes, we need to
        % be able to refer to vertex 1, 2 and 3.

        % Unfortunately, the only way I can think to do this is with a reverse
        % lookup loop, which isn't too efficient.  Hopefuly the small
        % searchlight patch will be fast enough, though.
        searchlightPatch_vis = [];
        for slVertex = searchlightPatchVs
            slVertex_i = find(indexMask.vertices == slVertex);
            searchlightPatch_vis = [searchlightPatch_vis, slVertex_i];
        end%for

        % Search through time
        window_i = 0;
        for window = slSpec.windowPositions'
            % thisWindow is the indices of timepoints in each window
            thisWindow = window(1):window(2);
            window_i = window_i + 1;

            searchlightPatchData = maskedMeshes(searchlightPatch_vis, thisWindow, :, :); % (vertices, time, condition, session)

            switch lower(userOptions.searchlightPatterns)
                case 'spatial'
                    % Spatial patterns: median over time window
                    searchlightPatchData = median(searchlightPatchData, 2); % (vertices, 1, conditions, sessions)
                    searchlightPatchData = squeeze(searchlightPatchData); % (vertices, conditions, sessions);
                case 'temporal'
                    % Temporal patterns: mean over vertices within searchlight
                    searchlightPatchData = mean(searchlightPatchData, 1); % (1, timePoints, conditions, sessions)
                    searchlightPatchData = squeeze(searchlightPatchData); % (timePionts, conditions, sessions)
                case 'spatiotemporal'
                    % Spatiotemporal patterns: all the data concatenated
                    searchlightPatchData = reshape(searchlightPatchData, [], size(searchlightPatchData, 3), size(searchlightPatchData, 4)); % (dataPoints, conditions, sessions)
            end%switch:userOptions.sensorSearchlightPatterns

            % Average RDMs over sessions
            searchlightRDM_tv = zeros(1, n_rdm_entries);
            for session = 1:nSessions
                sessionRDM = pdist(squeeze(searchlightPatchData(:,:,session))', userOptions.distance);
                searchlightRDM_tv = searchlightRDM_tv + sessionRDM;
            end%for:sessions
            searchlightRDM_tv = searchlightRDM_tv / nSessions;

            % Store results to be retured.
            searchlightRDMs(v_i, window_i, :) = searchlightRDM_tv(:);

        end%for:window

        % Indicate progress every once in a while...
        nVertsSearched = nVertsSearched + 1;
        if mod(nVertsSearched, 100) == 0
            prints('%d vertices searched', nVertsSearched);
        end%if

    end%for:v

    %% Saving RDM maps

    prints('Saving data RDMs to %s.', RDMsPath);
    save('-v7.3', RDMsPath, 'searchlightRDMs');
end%function
