userOptions = slOptions();


%% %%%%%%%%%%%%%%%%%%%
rsa.util.prints('Running toolbox for %s', userOptions.analysisName);
%%%%%%%%%%%%%%%%%%%%%%

dRDM = dynamic_bn_models();
n_lags = numel(dRDM);

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

averageRDMPaths = averageSlRDMs(RDMPaths, userOptions);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsa.util.prints( ...
   'Searchlight Brain Model comparison...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for subject_i = 1:nSubjects
   
   % Work on each hemisphere separately
   for chi = 'LR'
       
       rsa.util.prints();
       rsa.util.prints('Subject %d/%d (%sh)...', subject_i, n_subjects, chi);
       
       parfor lag_i = 1:n_lags
           
           this_lag_ms = (lag_i - 1) * MODEL_TIMESTEP_ms;
       
           rsa.util.prints('\tLag %d/%d...', lag_i, n_lags);

           % TODO: fix the use of sl metadatas and sl specs here
           [mapPaths(subject_i, lag_i).(chi)] = lag_fixed_searchlight_mapping_source( ...
               userOptions.subjectNames{subject_i}, ...
               chi, ...
               RDMPaths(subject_i).(chi), ...
               ...% Use the mask for this hemisphere only
               slMasks([slMasks.chi] == chi), ...
               dRDM(lag_i), this_lag_ms, ...
               STCMetadatas.(chi), ...
               userOptions ...
           );
   
       end
   end
end


%% Send an email

rsa.par.setupInternet();
rsa.par.setupEmail(userOptions.mailto);


%% Stop parallel toolbox

if userOptions.run_in_parallel;
    delete(p);
end
