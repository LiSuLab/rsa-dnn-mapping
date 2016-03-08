function [ averageRDMPaths ] = averageSlRDMs( RDMsPaths, slMasks, betaCorrs, userOptions )

    import rsa.*
    import rsa.util.*
    import rsa.par.*
    
    [nSubjects, nTimepoints] = size(RDMsPaths);
    
    % Paths
    file_i = 1;
    for chi = 'LR'
        for t = 1:nTimepoints
            averageRDMPaths(t).(chi) = fullfile( ...
                userOptions.rootPath, ...
                'RDMs', ...
                sprintf('average_t%d-%sh.mat', t, lower(chi)));
            promptOptions.checkFiles(file_i).address = averageRDMPaths(t).(chi);
            file_i = file_i + 1;
        end
    end
    
    promptOptions.functionCaller = 'averageSearchlightRDMs';
    promptOptions.defaultResponse = 'S';
    
    overwriteFlag = overwritePrompt(userOptions, promptOptions);
    
    if overwriteFlag

        for chi = 'LR'
            
            % We need these sizes to preallocate the matrices before
            % pushing off to the parfor.
            thisMask = slMasks([slMasks.chi] == chi);
            nVertices = numel(thisMask.vertices);
            nConditions = size(betaCorrs, 2);
            nEntries = numel(squareform(zeros(nConditions)));
                
            parfor t = 1:nTimepoints
                
                prints('Working on timepoint %d of %d...', t, nTimepoints);
            
                average_slRDMs = zeros(nVertices, nEntries);
                nan_counts = zeros(nVertices, nEntries);

                for subject_i = 1:nSubjects
    
                    this_subject_name = userOptions.subjectNames{subject_i};

                    prints('\tLoading searchlight RDMs for subject %s (%d/%d) %sh...', this_subject_name, subject_i, nSubjects, lower(chi));
                    this_subject_slRDMs = directLoad(RDMsPaths(subject_i, t).(chi));

                    % zero-out nans
                    nan_locations = isnan(this_subject_slRDMs);
                    this_subject_slRDMs(nan_locations) = 0;
                    % remember where the nans were
                    nan_counts = nan_counts + nan_locations;
                    % add to the average
                    average_slRDMs = average_slRDMs + this_subject_slRDMs;
                
                end%for:subject

                prints('\tAveraging RDMs at all vertices...');

                % replace nan counts by non-nan counts
                non_nan_counts = nSubjects - nan_counts;
                average_slRDMs = average_slRDMs ./ non_nan_counts;

                prints('\tSaving average searchlight RDMs for t=%d/%d to "%s"...', t, nTimepoints, averageRDMPaths(t).(chi));
                parsave(averageRDMPaths(t).(chi), average_slRDMs);
                
            end%for:t
        end%for:chi

    else
        prints('Average RDMs already calculated.  Skipping...');
    end

end

