import rsa.*
import rsa.util.*
import rsa.par.*
import rsa.meg.*

userOptions = bnMappingOptions();

% Here are some I made earlier
models = directLoad('/imaging/cw04/Neurolex/Lexpro/Analysis_DNN/Models/bn26_models.mat');

% The lag of the model timeline in miliseconds.
model_timeline_lag = ...
    ...% 100ms is the the offset for the alignment of the zero points.  We 
    ...% expect to see a fit for the models 100ms after the equivalent 
    ...% stimulus point in the brain data. This is consistent with the 
    ...% literature.
    100;


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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
   'Searchlight Brain Model comparison...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor subject_i = 1:nSubjects
   
   % Work on each hemisphere separately
   for chi = 'LR'
       rsa.util.prints('Working on subject %d, %sh side', subject_i, chi);

       % TODO: fix the use of sl metadatas and sl specs here
       [mapPaths(subject_i).(chi)] = rsa.meg.searchlight_dynamic_model_source( ...
           subject_i, ...
           chi, ...
           RDMPaths(subject_i).(chi), ...
           ...% Use the mask for this hemisphere only
           slMasks([slMasks.chi] == chi), ...
           model, ...
           [], ...
           adjacencyMatrices.L, ...
           STCMetadatas.(chi), ...
           userOptions ...
       );
   end
end


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
