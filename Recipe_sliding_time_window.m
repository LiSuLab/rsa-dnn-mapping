import rsa.*
import rsa.util.*
import rsa.par.*
import rsa.meg.*

userOptions = swMappingOptions();


%% %%%%%%%%%%%%%%%%%%%
prints('Running toolbox for %s', userOptions.analysisName);
%%%%%%%%%%%%%%%%%%%%%%

models = directLoad('/imaging/cw04/Neurolex/Lexpro/Analysis_DNN/Models/multilayer_model_RDMs_static_frame_5.mat');
n_models = numel(models);


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

swRDMPaths = slidingWindowRDMs_source(meshPaths, userOptions);

averageSWRDMPaths = averageSlidingWindowRDMs(swRDMPaths, userOptions);

for model_i = 1:n_models
    model = modelRDMs(model_i);
    
    prints('Sliding-window RSA for model "%s"...', model);
    
    % TODO: This shouldn't enforce left and right - it should be implicit from the mask name?
    ROI_sliding_TimeWindow(meshPaths, model, userOptions)
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
