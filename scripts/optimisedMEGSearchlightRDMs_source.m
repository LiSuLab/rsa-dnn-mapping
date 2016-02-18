function [RDMsPaths, slSTCMetadatas, slSpecs] = optimisedMEGSearchlightRDMs_source(meshPaths, slMasks, adjacencyMatrix, STCMetadatas, userOptions)

    import rsa.*
    import rsa.meg.*
    import rsa.rdm.*
    import rsa.util.*
    import rsa.par.*

    %% Constants

    nSubjects = numel(userOptions.subjectNames);

    [slSpecs, slSTCMetadatas] = getSearchlightSpec(STCMetadatas, userOptions);

    % The number of positions the sliding window will take.
    % There should be the same for both hemispheres
    nWindowPositions = size(slSpecs.L.windowPositions, 1);

    %% File paths

    RDMsDir = fullfile(userOptions.rootPath, 'RDMs');

    file_i = 1;
    for subject_i = 1:nSubjects
        thisSubjectName = userOptions.subjectNames{subject_i};
        for chi = 'LR'
            for t = 1:nWindowPositions
                RDMsPaths(subject_i, t).(chi) = fullfile( ...
                    RDMsDir, ...
                    sprintf('searchlightRDMs_%s_t%d-%sh.mat', thisSubjectName, t, lower(chi)));

                % We'll check all the files to be saved to see if they have already
                % been saved.
                promptOptions.checkFiles(file_i).address = RDMsPaths(subject_i, t).(chi);
                file_i = file_i + 1;
            end
        end
    end

    promptOptions.functionCaller = 'MEGSearchlightRDMs_source';
    promptOptions.defaultResponse = 'S';

    overwriteFlag = overwritePrompt(userOptions, promptOptions);


    %% Apply searchlight

    parfor subject_i = 1:nSubjects
        thisSubjectName = userOptions.subjectNames{subject_i};

        % Work on each hemisphere separately
        for chi = 'LR'

            prints('Shining RSA searchlight in the %sh source mesh of subject %d of %d (%s)...', lower(chi), subject_i, nSubjects, thisSubjectName);
            
            % If all files exist, we don't want to load a bunch of stuff if
            % we don't have to.
            all_files_exist = true;
            for t = 1:nWindowPositions
                if ~exist(RDMsPaths(subject_i, t).(chi), 'file')
                    all_files_exist = false;
                    break;
                end
            end
            if all_files_exist && ~overwriteFlag
                prints('All files exist for subject %s''s %sh. Skipping.', thisSubjectName, lower(chi));
                continue;
            end
            
            % Assuming we got this far there's something to do for this
            % subject's hemisphere.
            
            slSpec = slSpecs.(chi);
            indexMask = slMasks([slMasks.chi] == chi);
            meshPath = meshPaths(subject_i).(chi);

            maskedMeshes = directLoad(meshPath, 'sourceMeshes');

            [nVertices_masked, nTimePoints_data, nConditions, nSessions] = size(maskedMeshes);

            n_rdm_entries = numel(vectorizeRDM(zeros(nConditions)));

            for t = 1:nWindowPositions

                % We'll only do the searchlight if we haven't already done so,
                % unless we're being told to overwrite.
                if ~overwriteFlag && exist(RDMsPaths(subject_i, t).(chi), 'file')
                    prints('%s already exists; skipping...', RDMsPaths(subject_i, t).(chi));
                    continue;
                end

                prints('Working on timepont %d of %d...', t, nWindowPositions);

                RDMsPath = RDMsPaths(subject_i, t).(chi);

                %% map the volume with the searchlight

                % Preallocate looped matrices for speed
                searchlightRDMs = nan(nVertices_masked, n_rdm_entries);

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
                    window_limits = slSpec.windowPositions(t, :)';

                    % thisWindow is the indices of timepoints in each window
                    thisWindow = window_limits(1):window_limits(2);

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
                    searchlightRDMs(v_i, :) = searchlightRDM_tv(:);

                end%for:v

                %% Saving RDM maps

                prints('Saving data RDMs to %s.', RDMsPath);
                parsave(RDMsPath, searchlightRDMs);

            end%for:t

            %% Done
            prints('Done with subeject %d''s %sh side.', subject_i, lower(chi));

        end%for:chi
    end%for:subject

end%function
