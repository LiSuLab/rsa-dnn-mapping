function betas = lexproBetaCorrespondence()
%
%  betaCorrespondence.m is a simple function which should combine
%  three things: preBeta:	a string which is at the start of each file
%  containing a beta image, betas:	a struct indexed by (session,
%  condition) containing a sting unique to each beta image, postBeta:	a
%  string which is at the end of each file containing a beta image, not
%  containing the file .suffix
% 
%  use "[[subjectName]]" as a placeholder for the subject's name as found
%  in userOptions.subjectNames if necessary For example, in an experment
%  where the data from subject1 (subject1 name)  is saved in the format:
%  subject1Name_session1_condition1_experiment1.img and similarly for the
%  other conditions, one could use this function to define a general
%  mapping from experimental conditions to the path where the brain
%  responses are stored. If the paths are defined for a general subject,
%  the term [[subjectName]] would be iteratively replaced by the subject
%  names as defined by userOptions.subjectNames.
% 
%  note that this function could be replaced by an explicit mapping from
%  experimental conditions and sessions to data paths.
% 
%  Cai Wingfield 1-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

betaFilePath = '/imaging/at03/NKG_Data_Sets/LexproMEG/scripts/Stimuli-Lexpro-MEG-Single-col.txt';
betaNames = sort(textread(betaFilePath,'%s'));
preBeta = '[[subjectName]]-';
postBeta = '-[[LR]]h.stc';

nConditions = size(betaNames,1);

nTrials = size(betaNames,1)/nConditions;
%% betas
% all 400
for cond = 1:nConditions
    for trial = nTrials
        betas(trial,cond).identifier = [preBeta betaNames{(cond-1)*nTrials+trial} postBeta];     
    end
end

end%function
