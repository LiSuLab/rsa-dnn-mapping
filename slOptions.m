function userOptions = slOptions()

%%%%%%%%%%%%%%%%%%%%%
%% Project details %%
%%%%%%%%%%%%%%%%%%%%%

% This name identifies a collection of files which all belong to the same run of a project.
userOptions.analysisName = 'lexpro-bn-searchlight';

% This is the root directory of the project.
userOptions.rootPath = '/imaging/cw04/Neurolex/Lexpro/Analysis_DNN/CWD_win25_lateral';

% The path leading to where the scans are stored (not including subject-specific identifiers).
% "[[subjectName]]" should be used as a placeholder to denote an entry in userOptions.subjectNames
% "[[betaIdentifier]]" should be used as a placeholder to denote an output of betaCorrespondence.m if SPM is not being used; or an arbitrary filename if SPM is being used.
userOptions.betaPath = '/imaging/at03/NKG_Code_output/Version4_2/LexproMEG/3-single-trial-source-data/vert10242-smooth5-nodepth-eliFM-snr1-signed/[[betaIdentifier]]';

%%%%%%%%%%%%%%%%%%%
%% Email Options %%
%%%%%%%%%%%%%%%%%%%

% Set to true to be informed when a script finishes.
userOptions.recieveEmail = true;
% Put your address to make 
userOptions.mailto = 'cw417@cam.ac.uk';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parallel Computing toolbox %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To run parallel locally set run_in_parallel *and*
% run_in_parallel_in_cluster to false. 
% To run on CBU adaptive queue set run_in_parallel_in_cluster to true. 
% Do NOT set this true for fixed effect analysis (searchlight and
% permutation.
userOptions.run_in_parallel = true;
userOptions.run_in_parallel_in_cluster= true;
% Sometimes the performance will drop if there are large number of tiny
% jobs due to the communication and setting-up overhead.
userOptions.jobSize = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Adaptive Computing Cluster Queueing OPTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set true to delete jobs from the queue otherwise set to false.
% If you do not want to delete all jobs but only sepecific one do not set
% this variable to true.
userOptions.flush_Queue = true; 
% Used only when using CBU cluster.
% i.e. when run_in_parallel_in_cluster = true;
userOptions.wallTime = '24:00:00';
% Cluster machines requested.
userOptions.nodesReq = 6;
% Processors requested per processor machine.
userOptions.proPNode = 1;
% The product of nodesReq and proPNode should be greater or equal to the
% number of workers requested.
userOptions.nWorkers = 6;
% In gigabytes, to be distributed amongst all nodes.
userOptions.memReq = 270;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Modality-agnostic analysis options %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The path to a stereotypical mask data file is stored (not including subject-specific identifiers).
% "[[subjectName]]" should be used as a placeholder to denote an entry in userOptions.subjectNames
% "[[maskName]]" should be used as a placeholder to denote an entry in userOptions.maskNames
%userOptions.maskPath = '/imaging/fj01/cw04/New_Labels/[[maskName]]';
userOptions.maskPath = '/imaging/ef02/lexpro/subject/average/label/[[maskName]]';

% The list of mask filenames (minus .hdr extension) to be used.
% For MEG, names should be in pairs, such as maskName-lh,
% maskName-rh.
% Leave empty to do whole-brain analysis.
%
% For MEG sensor-level analysis, only the use of a single mask is
% supported.
userOptions.maskNames = { ...
    ...%'STG_STS_HG-lh', 'STG_STS_HG-rh', ...
    'lateral-lh', 'lateral-rh',...
};

% The type of pattern to look at.
% Options are:
%     Correlate over space ('spatial')
%     Correlate over time ('temporal')
%     Correlate over space and time ('spatiotemporal')
% For fMRI, the available options are 'spatial'.
% For MEG, the all options are available.
userOptions.searchlightPatterns = 'spatiotemporal';

%%%%%%%%%%%%%%%%%%%%%%%%
%% EXPERIMENTAL SETUP %%
%%%%%%%%%%%%%%%%%%%%%%%%

% The list of subjects to be included in the study.
userOptions.subjectNames = { ...
    'meg08_0320', ...      
    'meg08_0323', ...
    'meg08_0324', ...	
    'meg08_0327', ...
    'meg08_0348', ...		
    'meg08_0350', ... 
    'meg08_0363', ...
    'meg08_0366' ...
    'meg08_0371', ...
    'meg08_0372', ...    
    'meg08_0377' , ...
    'meg08_0380', ...	
    'meg08_0397', ...
    'meg08_0400', ...		
    'meg08_0401', ...
    'meg08_0402' ...
};% eg CBUXXXXX

% The default colour label for RDMs corresponding to RoI masks (as opposed to models).
userOptions.RoIColor = [0 0 1];
userOptions.ModelColor = [0 1 0];

%% %% %% %% %%
%%  MEG  %% Use these next four options if you're working in MEG:
%% %% %% %% %%

% The average surface files
userOptions.averageSurfaceFiles.L = '/imaging/ef02/lexpro/subject/average/surf/lh.inflated';
userOptions.averageSurfaceFiles.R = '/imaging/ef02/lexpro/subject/average/surf/rh.inflated';

% The width of the sliding window (ms)
userOptions.temporalSearchlightWidth = 25;

% The timestep for sliding window (ms)
userOptions.temporalSearchlightTimestep = 10;

% The overall window of interest for searchlight (ms)
userOptions.temporalSearchlightLimits = [0, 400];

% TODO: THis should n't be needed, bu tit is.
userOptions.maskTimeWindows = {[0,400], [0,400]};

% Temporal downsampling
% E.g., a value of 10 here means only taking each 10th point in time.
userOptions.temporalDownsampleRate = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MEG SOURCE-LEVEL ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The radius of the source-space searchlight (in mm)
userOptions.sourceSearchlightRadius = 20;

% Spatial downsampling.
% Set the target number of vertices per hemisphere.
userOptions.targetResolution = 5121;%10242;

% TODO: Explain this
% 5mm is the smallest distance between two adjacent vertex in 10242 resolution.
% 10mm is the smallest distance between two adjacent vertex in 2562 resolution.
userOptions.minDist = 5; %mm

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% First-order analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Text lables which may be attached to the conditions for MDS plots.
userOptions.conditionLabels = {};

% What colours should be given to the conditions?
userOptions.conditionColours = [repmat([1 0 0], 48,1); repmat([0 0 1], 44,1)];

% Which distance measure to use when calculating first-order RDMs.
userOptions.distance = 'Correlation';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Second-order analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Which similarity-measure is used for the second-order comparison.
userOptions.RDMCorrelationType = 'Spearman';

% How many permutations should be used to test the significance of the
% fits?  (10,000 highly recommended.)
userOptions.significanceTestPermutations = 10000;

% Bootstrap options
userOptions.nResamplings = 1000;
userOptions.resampleSubjects = true;
userOptions.resampleConditions = false;

userOptions.fisher = true;

% Group statistics options: random effect ('RFX') or fixed effect ('FFX')
% TODO: requires a lot more explanation
userOptions.groupStats = 'RFX';

% Clustering analysis: primary cluster-forming threshold in terms of the
% top x% for all vertexes across space and time.
% It means that primary threshold is set to p<x ,one tailed for positive
% clusters. This doesnot need to be adjusted according to the degree of
% freedom of data. If the requirement is to select top ONLY x% of all data, set
% this value to nan.
% If not using tmaps for RFX, this value only will define the primary threshold.
% DoF does not help compute threshold for that case.
userOptions.primaryThreshold = 0.05;

% Should RDMs' entries be rank transformed into [0,1] before they're displayed?
userOptions.rankTransform = true;

%%%%%%%%%%%%%%%%%%%%
%% Figure options %%
%%%%%%%%%%%%%%%%%%%%

% Should rubber bands be shown on the MDS plot?
userOptions.rubberbands = false;

% What criterion shoud be minimised in MDS display?
userOptions.criterion = 'metricstress';

% What is the colourscheme for the RDMs?
userOptions.colourScheme = hot;

% How should figures be outputted?
userOptions.displayFigures = true;
userOptions.saveFiguresPDF = false;
userOptions.saveFiguresFig = false;
userOptions.saveFiguresPS = false;
% Which dots per inch resolution do we output?
userOptions.dpi = 300;
% Remove whitespace from PDF/PS files?
% Bad if you just want to print the figures since they'll
% no longer be A4 size, good if you want to put the figure
% in a manuscript or presentation.
userOptions.tightInset = false;

%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interaction options %%
%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Need a better solution to this. Enum?
% This can be used to force a reply to each propt about overwriting files.
% It can be useful when running unsupervised, or in parallel, where the
% prompt may not even be seen.
% Set to 'a', 'r', 's'; or use '' to not force a reply.
userOptions.forcePromptReply = 's';

end%function
