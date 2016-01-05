userOptions = swMappingOptions();


%% %%%%%%%%%%%%%%%%%%%
rsa.util.prints('Running toolbox for %s', userOptions.analysisName);
%%%%%%%%%%%%%%%%%%%%%%

modelRDMs = rsa.util.directLoad('/imaging/cw04/Neurolex/Lexpro/Analysis_DNN/Models/multilayer_model_RDMs_static_frame_5.mat');
n_models = numel(modelRDMs);


%% %%%%%%%%%%%%%%%%%%%
rsa.util.prints('Preparing masks...');
%%%%%%%%%%%%%%%%%%%%%%

% TODO: Don't enforce use of both hemispheres

if ~isempty(userOptions.maskNames)
    slMasks = rsa.meg.MEGMaskPreparation_source(userOptions);
    % For this searchlight analysis, we combine all masks into one
    slMasks = rsa.meg.combineVertexMasks_source(slMasks, 'combined_mask', userOptions);  
else
    slMasks = rsa.meg.allBrainMask(userOptions);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints('Starting parallel toolbox...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if userOptions.flush_Queue
    rsa.par.flushQ();
end

if userOptions.run_in_parallel
    p = rsa.par.initialise_CBU_Queue(userOptions);
end


%% %%%%%%%%%%%%%%%%%%
rsa.util.prints('Loading brain data...');
%%%%%%%%%%%%%%%%%%%%%

% TODO: Bring the parfor loop outside of this

[meshPaths, STCMetadatas] = rsa.meg.MEGDataPreparation_source( ...
    lexproBetaCorrespondence(), ...
    userOptions, ...
    'mask', slMasks);

swRDMsPaths = rsa.meg.MEGSlidingWindowRDMs_source(meshPaths, STCMetadatas, userOptions);

aswRDMsPaths = rsa.meg.averageSlidingWindowRDMs(swRDMPaths, userOptions);

for model_i = 1:n_models
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
