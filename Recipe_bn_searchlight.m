userOptions = slOptions_hidden_layer();


%% %%%%%%%%%%%%%%%%%%%
rsa.util.prints('Running toolbox for %s', userOptions.analysisName);
%%%%%%%%%%%%%%%%%%%%%%

dynamic_model_RDM = dynamic_hidden_layer_models('correlation');
n_lags = numel(dynamic_model_RDM);

MODEL_TIMESTEP_ms = 10;


%% %%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'Preparing masks...');
%%%%%%%%%%%%%%%%%%%%%%

usingMasks = ~isempty(userOptions.maskNames);
if usingMasks
    slMasks = rsa.meg.MEGMaskPreparation_source(userOptions);
    % For this searchlight analysis, we combine all masks into one
    slMasks = rsa.meg.combineVertexMasks_source(slMasks, 'combined_mask', userOptions);  
else
    slMasks = rsa.meg.allBrainMask(userOptions);
end


%% Compute some constats
nSubjects = numel(userOptions.subjectNames);
adjacencyMatrix = rsa.meg.calculateMeshAdjacency(userOptions.targetResolution, userOptions.sourceSearchlightRadius, userOptions);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'Starting parallel toolbox...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if userOptions.flush_Queue
    rsa.par.flushQ();
end

if userOptions.run_in_parallel
    p = rsa.par.initialise_CBU_Queue(userOptions);
end


%% %%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'Loading brain data...');
%%%%%%%%%%%%%%%%%%%%%

[meshPaths, STCMetadatas] = rsa.meg.MEGDataPreparation_source( ...
    lexproBetaCorrespondence(), ...
    userOptions, ...
    'mask', slMasks);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
    'Searchlight Brain RDM Calculation...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[RDMsPaths, slSTCMetadatas, slSpecs] = optimisedMEGSearchlightRDMs_source( ...
    meshPaths, ...
    slMasks, ...
    adjacencyMatrix, ...
    STCMetadatas, ...
    userOptions);


%% %%%%%
rsa.util.prints( ...
    'Averaging searchlight RDMs...');
%%%%%%%%

averageRDMPaths = averageSlRDMs( ...
    RDMsPaths, ...
    slMasks, ...
    lexproBetaCorrespondence(), ...
    userOptions);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
   'Searchlight Brain Model comparison...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_subjects = numel(userOptions.subjectNames);
   
% Work on each hemisphere separately
for chi = 'LR'
    
    for subject_i = 1:n_subjects
        
        subjectName = userOptions.subjectNames{subject_i};

        for lag_i = 1:n_lags

           this_lag_ms = (lag_i - 1) * MODEL_TIMESTEP_ms;

           rsa.util.prints('%sH: Subj %d/%d. Lag %d/%d...', chi, subject_i, n_subjects, lag_i, n_lags);

           % TODO: fix the use of sl metadatas and sl specs here
           [mapPaths(subject_i, lag_i).(chi)] = lag_fixed_searchlight_mapping_source( ...
               chi, ...
               ...% can't deref .(chi) here as it auto-expands the argument
               RDMsPaths(subject_i, :), ...
               sprintf('%s_subj_%s', userOptions.analysisName, subjectName), ...
               ...% Use the mask for this hemisphere only
               slMasks([slMasks.chi] == chi), ...
               dynamic_model_RDM(lag_i), this_lag_ms, ...
               ...% Apparently we don't deref .(chi) here either
               STCMetadatas, ...
               userOptions);

        end
    end
end

% Work on subject average
for chi = 'LR'
    for lag_i = 1:n_lags
        this_lag_ms = (lag_i - 1) * MODEL_TIMESTEP_ms;
        [gaMapPaths(lag_i).(chi)] = lag_fixed_searchlight_mapping_source( ...
            chi, ...
            ...% can't deref .(chi) here as it auto-expands the argument
            averageRDMPaths, ...
            sprintf('%s_groupAverageRDM', userOptions.analysisName), ...
            ...% Use the mask for this hemisphere only
            slMasks([slMasks.chi] == chi), ...
            dynamic_model_RDM(lag_i), this_lag_ms, ...
            ...% Apparently we don't deref .(chi) here either
            STCMetadatas, ...
            userOptions);
    end
end%for


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints(...
    'Averaging STC files...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for subject_i = 1:n_subjects
    
    subjectName = userOptions.subjectNames{subject_i};
   
    integratedMapPaths(subject_i) = lag_integrate_stc_files(...
        mapPaths(subject_i, :), ...
        sprintf('%s_subj_%s_lagfixed', userOptions.analysisName, subjectName), ...
        userOptions, ...
        true);

end

gaIntegratedMapPaths = lag_integrate_stc_files(...
    gaMapPaths, ...
    sprintf('%s_groupAverageRDM_lagfixed', userOptions.analysisName), ...
    userOptions, ...
    true);

averageMapPaths = average_stc_files( ...
    integratedMapPaths, ...
    sprintf('%s_subave_lagfixed', userOptions.analysisName), ...
    userOptions);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints(...
    'RFX stats...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cluster_forming_threshold = 0.01;
fdr_threshold = 0.05;
n_permutations = 1000;

rsa.util.prints('Simulating statistial-maps...');

[observed_map_paths, corrected_ps] = rfx_cluster( ...
    integratedMapPaths, ...
    n_permutations, ...
    ...% statistic type
    'r', ...
    cluster_forming_threshold, ...
    fdr_threshold, ...
    userOptions);



%% Send an email

rsa.par.setupInternet();
rsa.par.setupEmail(userOptions.mailto);


%% Stop parallel toolbox

if userOptions.run_in_parallel;
    delete(p);
end
