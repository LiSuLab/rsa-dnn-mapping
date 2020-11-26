% Expects phonetic_model_RDMs to be in the workspace.
% Run CalculateAndShowPhoneRDMs.m first.

FEATURES = phonetic_feature_matrix_GMM();

feature_names = fields(FEATURES);

n_features = numel(feature_names);

[n_timepoints, n_models] = size(phonetic_model_RDMs);

skip_frames = 4;

% stack RDMs ahead of time
clear RDM_data;
RDM_data(n_timepoints) = struct();
for t = 1:n_timepoints
    RDM_data(t).data = [];
    for m = 1:n_models
        RDM_data(t).data = [RDM_data(t).data; phonetic_model_RDMs(t, m).RDM'];
    end
end

n_timepoints = n_timepoints - skip_frames;

db_is = nan(n_features, n_timepoints);

parfor t = 1:n_timepoints

    rsa.util.prints('Timepoint %d of %d...', t, n_timepoints);

    for feature_i = 1:n_features

        feature_name = feature_names{feature_i};

        rsa.util.prints('\tFeature %d of %d...', feature_i, n_features);

        % For each feature we'll consider the between-feat/nonfeat distances
        % and the within feat/nonfeat distances.

        % 0s and 1s
        feature_profile = FEATURES.(feature_name);
        
        % 1s and 2s
        feature_class_is = feature_profile + 1;
        
        [db_is(feature_i, t), max_ignore] = db_index(RDM_data(t).data, feature_class_is);

    end%for:feature_i

end%for:t

rsa.util.prints('Done!');
