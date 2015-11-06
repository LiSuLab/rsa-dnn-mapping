% TODO: Documentation
%
% Cai Wingfield 2010-05, 2010-08, 2015-03--2015-06
% update by Li Su 3-2012, 11-2012
% updated Fawad 12-2013, 02-2014, 10-2014

import rsa.*
import rsa.util.*
import rsa.par.*
import rsa.meg.*

userOptions = phoneticMappingOptions();

prints('Starting RSA analysis "%s".', userOptions.analysisName);


%% %%%%%%%%%%%%%%%%%%%%%%%%
prints('Preparing model RDMs...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%models = constructModelRDMs(userOptions);
% Here are some I made earlier
models = directLoad('/imaging/cw04/Neurolex/Lexpro/Analysis_Phonetic_mapping/Model_HTK_dnn/triphone-likelihood-RDMs.mat');

% Trim the unusable frames from the beginning of the model timeline.
trim_frames = 4;
models = models(trim_frames+1:end, :);

% The lag of the model timeline in miliseconds.
model_timeline_lag = ...
    ...% 100ms is the the offset for the alignment of the zero points.  We 
    ...% expect to see a fit for the models 100ms after the equivalent 
    ...% stimulus point in the brain data. This is consistent with the 
    ...% literature.
    100 ...
    ...% We further increase the lag to account for the fact that we're
    ...% trimming frames from the start of the model timeline, and each
    ...% frame is 10ms long.
    + (trim_frames * 10);


%% %%%%%%%%%%%%%%%%%%%
prints('Preparing masks...');
%%%%%%%%%%%%%%%%%%%%%%

% TODO: Don't enforce use of both hemispheres

usingMasks = ~isempty(userOptions.maskNames);
if usingMasks
    slMasks = MEGMaskPreparation_source(userOptions);
    % For this searchlight analysis, we combine all masks into one
    slMasks = combineVertexMasks_source(slMasks, 'combined_mask', userOptions);  
else
    slMasks = allBrainMask(userOptions);
end

adjacencyMatrices = calculateMeshAdjacency(userOptions.targetResolution, userOptions.sourceSearchlightRadius, userOptions, 'hemis', 'LR');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
prints('Starting parallel toolbox...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if userOptions.flush_Queue
    flushQ();
end

if userOptions.run_in_parallel
    p = initialise_CBU_Queue(userOptions);
end


%% %%%%%%%%%%%%%%%%%%
prints('Loading brain data...');
%%%%%%%%%%%%%%%%%%%%%

[meshPaths, STCMetadatas] = MEGDataPreparation_source( ...
    lexproBetaCorrespondence(), ...
    userOptions, ...
    'mask', slMasks);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prints('Searchlight Brain RDM Calculation...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[RDMsPaths, slSTCMetadatas] = MEGSearchlightRDMs_source( ...
    meshPaths, ...
    slMasks, ...
    ...% Assume that both hemis' adjacency matrices are the same so only use one.
    adjacencyMatrices.L, ...
    STCMetadatas, ...
    userOptions);


%% %%%%%
prints('Averaging searchlight RDMs...');
%%%%%%%%

averageRDMPaths = averageSearchlightRDMs(RDMsPaths, userOptions);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prints('GLM-fitting models to searchlight RDMs...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[glm_paths, lagSTCMetadatas] = searchlight_dynamicGLM_source( ...
    averageRDMPaths, ...
    models, ...
    slSTCMetadatas, ...
    userOptions, ...
    'lag', model_timeline_lag);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prints('Thresholding GLM values...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h0_paths] = searchlight_GLM_permutation_source( ...
    averageRDMPaths, ...
    models, ...
    slSTCMetadatas, ...
    lagSTCMetadatas, ...
    ...% number of permutation batches
    ...%TODO make this into optional argument
    30, ... % 30
    userOptions);
    
FEATURES = phonetic_feature_matrix_dnn();

[feature_paths, feature_paths_ea] = label_with_features( ...
    glm_paths, ...
    FEATURES, ...
    lagSTCMetadatas, ...
    true, ...
    userOptions);

% Compute thresholds for different levels of significance.
% Index equals number of stars.
separate_fit_thresholds = null_distribution_of_sums( ...
    h0_paths, ...
    FEATURES, ...
    'threshold', 0.05);

% Print out the thresholded values
rsa.util.display_singleton_struct(separate_fit_thresholds);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prints('Anlysing feature patches...');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Check to see if the thresholded features in the paper are being extracted
% % properly.

% feature_threshold_levels = custom_feature_thresholds_dnn();
% 
% [feature_thresholded_paths, threshold_vis] = threshold_feature_maps( ...
%     feature_paths_ea, ...
%     feature_threshold_levels, ...
%     separate_fit_thresholds, ...
%     userOptions);
%     
% % TODO: get this path right, maybe include the to-be-loaded file it in the repo?
% tonotopy_paths = directLoad('/imaging/cw04/Neurolex/Lexpro/Analysis_Phonetic_mapping/analysis-phonetc-mapping-scripts/tonotopy_paths.mat');
% 
% [beta_responses, beta_responses_positive, beta_sds] = group_tonotopy_betas_inside_patches( ...
%     threshold_vis, ...
%     tonotopy_paths);
% 
% save_beta_figues( ...
%     beta_responses, ...
%     userOptions);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prints('Cleaning up...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Close the parpool
if userOptions.run_in_parallel
    delete(p);
end

% Sending an email
if userOptions.recieveEmail
    setupInternet();
    setupEmail(userOptions.mailto);
end

prints( ...
    'RSA COMPLETE!');
