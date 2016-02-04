function [ averageRDMPaths ] = averageSlRDMs( RDMPaths, userOptions )

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

            for subject_i = 1:nSubjects

                this_subject_name = userOptions.subjectNames{subject_i};

                prints('Loading searchlight RDMs for subject %s (%d/%d) %sh...', this_subject_name, subject_i, nSubjects, lower(chi));
                this_subject_slRDMs = directLoad(RDMPaths(subject_i).(chi), 'searchlightRDMs');

                % For the first subject, we initialise the average and the
                % nan-counter with some sizes.
                if subject_i == 1
                    [nVertices, nTimepoints, nEntries] = size(this_subject_slRDMs);
                    average_slRDMs = nan(nVertices, nTimepoints, nEntries);
                    nan_counts = zeros(1:nVertices, 1:nTimepoints, nEntries);
                end

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

