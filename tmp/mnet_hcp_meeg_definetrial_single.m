function [D,fieldtripData] = mnet_hcp_meeg_definetrial_single(config, D, fieldtripData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEEG Redefine Trials According to Band Power.                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.06.16 02:10 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input
% These five fields must be specified
if ~isfield(config, 'subjectID'),error('Invlid Subject ID!! There is no subjectID');end
if (~isfield(config, 'subjectPath')||~isdir(config.subjectPath));error('Invalid Subject Path!! There is no such folder');end
if (~isfield(config, 'savePath')||~isdir(config.savePath)),error('Invalid Save Path!! There is no such folder');end
if (~isfield(config, 'sessionID')),error('Invalid Session ID!! There is no session ID field!!');end
%% Set Path
rMEGPath = fullfile(config.subjectPath,[config.subjectID '_MEG_Restin_preproc'],config.subjectID,'MEG','Restin','rmegpreproc');
rMEGRawData=[config.subjectID '_MEG_' num2str(config.sessionID) '-Restin_rmegpreproc'];
anaPath = fullfile(config.subjectPath, [config.subjectID '_MEG_anatomy'],config.subjectID, 'MEG','anatomy');
headmodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_headmodel']);
sourcemodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_sourcemodel_2d']);
%% Load MEEG Datas
use_my_library('ft',0);
use_my_library('spm',1);
if nargin==1 || isempty(D)||isempty(fieldtripData)
    try
        D = spm_eeg_load(fullfile(config.savePath,['affdspm_', rMEGRawData]));
        fieldtripData = load(fullfile(rMEGPath,rMEGRawData));
    catch
        error(['There is no preprocessed file!! filename : ' fullfile(config.savePath,['affdspm_', rMEGRawData])]);
    end
end
% if length(D.condlist)>1
%     condpart = strsplit(D.condlist{1}, ' ');
%     if strcmp(lower(condpart{1}),config.band)
%         return; % If already calculated, then return immediately.
%     end
% end
%% Redefine trials
conditionlabels = cell(1,D.ntrials);
for idx = 1:D.ntrials
    conditionlabels{idx} = num2str(idx);
end
D = conditions(D,':',conditionlabels);
save(fullfile(config.savePath,D.fname),'D');
