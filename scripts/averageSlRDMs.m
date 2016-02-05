function [ averageRDMPaths ] = averageSlRDMs( RDMPaths, slSpecs, slMasks, betaCorrs, userOptions )

    import rsa.*
    import rsa.util.*
    
    % Paths
    file_i = 1;
    for chi = 'LR'
        averageRDMPaths.(chi) = fullfile(userOptions.rootPath, 'RDMs', ['average_', lower(chi), 'h.mat']);
        promptOptions.checkFiles(file_i).address = averageRDMPaths.(chi);
        file_i = file_i + 1;
    end
    
    promptOptions.functionCaller = 'averageSearchlightRDMs';
    promptOptions.defaultResponse = 'S';
    
    overwriteFlag = overwritePrompt(userOptions, promptOptions);
    
    if overwriteFlag

        nSubjects = numel(userOptions.subjectNames);
        for chi = 'LR'
            
            % We need these sizes to preallocate the matrices before
            % pushing off to the parfor.
            thisMask = slMasks([slMasks.chi] == chi);
            nVertices = numel(thisMask.vertices);
            nTimepoints = size(slSpecs.(chi).windowPositions, 1);
            nConditions = size(betaCorrs, 2);
            nEntries = numel(squareform(zeros(nConditions)));
            
            average_slRDMs = nan(nVertices, nTimepoints, nEntries);
            nan_counts = zeros(nVertices, nTimepoints, nEntries);

            parfor subject_i = 1:nSubjects

                this_subject_name = userOptions.subjectNames{subject_i};

                prints('Loading searchlight RDMs for subject %s (%d/%d) %sh...', this_subject_name, subject_i, nSubjects, lower(chi));
                this_subject_slRDMs = directLoad(RDMPaths(subject_i).(chi), 'searchlightRDMs');

                prints('Adding RDMs at all vertices and timepoints...');
                
                % zero-out nans
                nan_locations = isnan(this_subject_slRDMs);
                this_subject_slRDMs(nan_locations) = 0;
                % remember where the nans were
                nan_counts = nan_counts + nan_locations;
                % add to the average
                average_slRDMs = average_slRDMs + this_subject_slRDMs;

            end%for:subject

            prints('Averaging RDMs at all vertices...');

            % replace nan counts by non-nan counts
            non_nan_counts = nSubjects - nan_counts;
            average_slRDMs = average_slRDMs ./ non_nan_counts;

            prints('Saving average searchlight RDMs to "%s"...', averageRDMPaths.(chi));
            save('-v7.3', averageRDMPaths.(chi), 'average_slRDMs');

        end%for:chi

    else
        prints('Average RDMs already calculated.  Skipping...');
    end

end

