userOptions = swMappingOptions();


%% %%%%%%%%%%%%%%%%%%%
prints('Running toolbox for %s', userOptions.analysisName);
%%%%%%%%%%%%%%%%%%%%%%

modelRDMs = directLoad('/imaging/cw04/Neurolex/Lexpro/Analysis_DNN/Models/multilayer_model_RDMs_static_frame_5.mat');
n_models = numel(modelRDMs);


%% %%%%%%%%%%%%%%%%%%%
prints('Preparing masks...');
%%%%%%%%%%%%%%%%%%%%%%

% TODO: Don't enforce use of both hemispheres

if ~isempty(userOptions.maskNames)
    slMasks = MEGMaskPreparation_source(userOptions);
    % For this searchlight analysis, we combine all masks into one
    slMasks = combineVertexMasks_source(slMasks, 'combined_mask', userOptions);  
else
    slMasks = allBrainMask(userOptions);
end

% TODO: Only need one mesh adjacency here - they're the same for both
% hemispheres.

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

% TODO: Bring the parfor loop outside of this

[meshPaths, STCMetadatas] = MEGDataPreparation_source( ...
    lexproBetaCorrespondence(), ...
    userOptions, ...
    'mask', slMasks);

swRDMsPaths = MEGSlidingWindowRDMs_source(meshPaths, STCMetadatas, userOptions);

aswRDMsPaths = averageSlidingWindowRDMs(swRDMPaths, userOptions);

parfor model_i = 1:n_models
    modelRDM = modelRDMs(model_i);
    
    prints('Sliding-window RSA for model "%s"...', modelRDM);
    
    % TODO: This shouldn't enforce left and right - it should be implicit from the mask name?
    output_Rs(model_i, :) = rsa.meg.sliding_time_window_source(aswRDMsPaths, modelRDM, userOptions);
end


%% Send an email

if userOptions.recieveEmail
    rsa.par.setupInternet();
    rsa.par.setupEmail(userOptions.mailto);
end


%% Stop parallel toolbox

if userOptions.run_in_parallel;
    delete(p);
end
