%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCP rMEG Intra-band Contrast spDCM (ERP model) Batch                    %
%                                                                         %
% HCP Resting State MEG Dataset has been cut into 2 second trials         %
% This pipeline split trials into two sets according to bandpower of      %
% interest. Then by using parametric empirical Baysian Methods, summarize %
% group-level parameter estimation.                                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.06.16 01:40 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Subject from HCP Dataset
clear;
addpath('/home/kuro/Codes/M_code/my_code'); % This is the directory where the ralated codes are located
use_my_library('spm',1);
mPath='/projects1/HCPMEG'; % This is the directory of root of HCP dataset.
workingDir = '/projects1/MEGResult'; % This is the directory that all result files will be placed.
megDirInfo = dir(mPath);
subject = [];
cd(mPath);
for idx=1:length(megDirInfo)
    if ~isdir(megDirInfo(idx).name)
        continue;
    end
    if ~isempty(strfind(megDirInfo(idx).name, 'MEG')) % if subject folder is '[HCP data path]/[subject id]_MEG'
        subject(end+1).id = strtok(megDirInfo(idx).name,'_');
        subject(end).dir = megDirInfo(idx).name;
    elseif ~isempty(str2num(megDirInfo(idx).name)) % if subject folder is '[HCP data path]/[subject id]'
        subject(end+1).id = megDirInfo(idx).name;
        subject(end).dir = megDirInfo(idx).name;
    end
end
clear megDirInfo;
for idx =1:length(subject)
    subjDirInfo = dir(fullfile(mPath,subject(idx).dir));
    subject(idx).anatomy = [];
    subject(idx).Restin = [];
    subject(idx).Wrkmem = [];
    subject(idx).StoryM = [];
    subject(idx).Motort = [];
    for j =1:length(subjDirInfo)
       if ~isempty(strfind(subjDirInfo(j).name, 'anatomy'))
           subject(idx).anatomy = 1; 
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Restin_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(idx).dir,subjDirInfo(j).name,subject(idx).id, 'MEG/Restin/rmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(idx).Restin(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(idx).Restin =unique(subject(idx).Restin); 
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Wrkmem_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(idx).dir,subjDirInfo(j).name,subject(idx).id, 'MEG/Wrkmem/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(idx).Wrkmem(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(idx).Wrkmem =unique(subject(idx).Wrkmem);
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'StoryM_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(idx).dir,subjDirInfo(j).name,subject(idx).id, 'MEG/StoryM/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(idx).StoryM(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(idx).StoryM =unique(subject(idx).StoryM);
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Motort_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(idx).dir,subjDirInfo(j).name,subject(idx).id, 'MEG/Motort/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(idx).Motort(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(idx).Motort =unique(subject(idx).Motort);
       end
    end
end
clear tempDirInfo subjDirInfo subDirInfo sessionNum i j k;

cd(workingDir);
restWorkingDir               = fullfile(workingDir, 'Restin');
if ~isdir(restWorkingDir),mkdir(restWorkingDir);end
dataset = cell(length(subject), 3, 2); % # of subjects x # of sessions x 2 datatype. One for spm MEEG data, another for fieldtrip data.
cfg.savePath                 = restWorkingDir;
%% Load and Do Preprocess
cfg.freqband                 = [0.1,40];
cfg.downsample               = 200;

for idx = 1:length(subject)
    cfg.subjectID            = subject(idx).id;
    cfg.subjectPath          = fullfile(mPath, subject(idx).dir);
    for sess = 1:length(subject(idx).Restin)
        cfg.sessionID        = subject(idx).Restin(sess);
        [dataset{idx,sess,1},dataset{idx,sess,2}] = mnet_hcp_meeg_preprocess_spm(cfg); % load and preprocess data.
    end
end
% You should check by these two fuctions
% D = spm_eeg_load('filename');
% spm_eeg_review(D);
%% Set Forward Model
for idx = 1:length(subject)
    cfg.subjectID            = subject(idx).id;
    cfg.subjectPath          = fullfile(mPath, subject(idx).dir);
    for sess = 1:length(subject(idx).Restin)
        cfg.sessionID        = subject(idx).Restin(sess);
        [dataset{idx,sess,1},dataset{idx,sess,2}] = mnet_hcp_meeg_forward(cfg,dataset{idx,sess,1},dataset{idx,sess,2});
    end
end
%% Redefine Trials As Single Trials
for idx = 1:length(subject)
    cfg.subjectID            = subject(idx).id;
    cfg.subjectPath          = fullfile(mPath, subject(idx).dir);
    for sess = 1:length(subject(idx).Restin)
        cfg.sessionID        = subject(idx).Restin(sess);
        dataset{idx,sess,1}=mnet_hcp_meeg_definetrial_single(cfg,dataset{idx,sess,1},dataset{idx,sess,2});
    end
end
%% Run spDCM
cfg.TimeWin                  = [0 2000];
cfg.FreqWin                  = [1,40];
cfg.analysis                 = 'CSD';
cfg.model                    = 'ERP';
cfg.spatial                  = 'ECD';
% cfg.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27] [-28; -96; -6] [28;-96;-6]];
% cfg.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC', 'lV1','rV1'};
cfg.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27]];
cfg.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC'};
cfg.chanMode                 = 8;

cfg.A{1}                     = [0 0 0 0; 0 0 0 0; 1 1 0 0 ; 1 1 1 0];      % Forward
cfg.A{2}                     = cfg.A{1}';                              % Backward
cfg.A{3}                     = [0 1 0 0 ;1 0 0 0; 0 0 0 0; 0 0 0 0];       % Modulatory

% Full Model
% cfg.A{1}                     = ones(4,4);                                    % Forward
% cfg.A{2}                     = ones(4,4);                                    % Backward
% cfg.A{3}                     = ones(4,4);                                    % Modulatory

cfg.B{1}                     = zeros(4,4);

cfg.effect                   = zeros(1,0); % Between trial effect (Design matrix)
cfg.effectName               = {};

if ~isdir(fullfile(cfg.savePath,'DCMs'))
    mkdir(fullfile(cfg.savePath,'DCMs'));
end

for idx = 1:length(subject)
    cfg.subjectID            = subject(idx).id;
    cfg.subjectPath          = fullfile(mPath, subject(idx).dir);
    for sess = 1:length(subject(idx).Restin)
        cfg.sessionID        = subject(idx).Restin(sess);
        for conditions = 1:dataset{idx,sess,1}.ntrials
            cfg.trials           = [conditions]; % D.condlist{1}, LowAlpha
            [D, DCM]   = mnet_hcp_meeg_rest_spDCM_erpmodel(cfg,dataset{idx,sess,1});
            save(fullfile(cfg.savePath,'DCMs',['DCM_' cfg.subjectID '_sess' num2str(cfg.sessionID) '_trial' num2str(conditions)]),'DCM');
            close all;
        end
    end
end

%% PEB
% % load DCMs
% dirInfo = dir(fullfile(restWorkingDir,'DCMs','Theta'));
% fileRemoveIdx = zeros(0,1);
% for idx =1:length(dirInfo)
%     fileRemoveIdx(end+1) = isempty(strfind(dirInfo(idx).name, 'DCM'))|| isempty(strfind(dirInfo(idx).name, [upper(cfg.band(1)) cfg.band(2:end)]));
% end
% dirInfo(find(fileRemoveIdx)) = []; % get filenames contain string 'DCM'.
% clear fileRemoveIdx i
% 
% Descriptions = [];                                                         
% DCMs = cell(0,2);                                                          
% 
% idx = 1;
% while idx < length(dirInfo)
%     fileNameSplit1 = strsplit(dirInfo(idx).name,'_');
%     fileNameSplit2 = strsplit(dirInfo(idx+1).name,'_');
%     if strcmp(fileNameSplit1{2},fileNameSplit2{2}) && strcmp(fileNameSplit1{3},fileNameSplit2{3})
%         Descriptions(end+1).subjectId = fileNameSplit1{2};
%         Descriptions(end).sessionId = fileNameSplit1{3};
%         load(fullfile(fullfile(restWorkingDir,'DCMs','Theta'),dirInfo(idx).name));
%         load(fullfile(fullfile(restWorkingDir,'DCMs','Theta'),dirInfo(idx+1).name));
%         DCMs(end+1,:) = {DCM_High,DCM_Low};
%     else
%         idx=idx+1;
%         continue;
%     end
%     idx = idx+2;
% end
% % Sort DCMs By sessionId
% Descriptions_new= [];   
% DCMs_old = DCMs;
% DCMs = cell(0,6);
% idx = 1;
% while idx < length(Descriptions)
%     if~(strcmp(Descriptions(idx).subjectId, Descriptions(idx+1).subjectId)&&...
%             strcmp(Descriptions(idx+1).subjectId,Descriptions(idx+2).subjectId))
%         idx=idx+1;
%         continue;
%     end
%     Descriptions_new(end+1).subjectId =Descriptions(idx).subjectId;
%     DCMs(end+1,:) = {DCMs_old{idx,1},DCMs_old{idx+1,1},DCMs_old{idx+2,1},DCMs_old{idx,2},DCMs_old{idx+1,2},DCMs_old{idx+2,2}};
%     idx = idx+3;
% end
% Descriptions = Descriptions_new;
% clear DCMs_old DCM_High DCM_Low Descriptions_new i dirInfo fileNameSplit1 fileNameSplit2
% 
% PEBi = {}; M1.X = ones(6,2); M1.X(4:6,2) = -1;
% for s = 1:size(DCMs, 1)
%     dcm = {};
%     for j=1:6, dcm{j,1}=DCMs{s,j};end
%     PEBi{s,1} = spm_dcm_peb(dcm,M1);
% end
% use_my_library('mnet',1);
% PEBa = spm_dcm_peb_peb(PEBi);
% mnet_dcm_peb_display(PEBa);
% use_my_library('mnet',0);
% %%
% freqbands = {'theta', 'delta', 'beta','gamma1','gamma2','alpha'};
% for j = 4:6
%     %% Redefine Trials According to Band Power
%     try
%         cfg                      = rmfield(cfg, 'freqband');
%         cfg                      = rmfield(cfg, 'downsample');
%     end
%     % set band : 'delta' (2-4Hz), 'theta' (5-7Hz), 'alpha' (8-12Hz), 
%     %    'beta' (15-30Hz), 'gamma1' (30-60Hz), 'gamma2(60-90Hz)
%     % Trials defined by median power of median frequency of each band. 
%     %    ex. for alpha, calculate 10Hz(median frequency of alpha) power for 
%     %        each trials and separate trials two group with larger than median
%     %        trial power and smaller than median trial power.
%     cfg.band                     = freqbands{j};
%     cfg.showPowerRatio           = 0;
%     for idx = 1:length(subject)
%         cfg.subjectID            = subject(idx).id;
%         cfg.subjectPath          = fullfile(mPath, subject(idx).dir);
%         for sess = subject(idx).Restin
%             cfg.sessionID        = sess;
%             mnet_hcp_meeg_definetrial_bandpower(cfg);
%         end
%     end
%     %% Set Forward Model
%     for idx = 1:length(subject)
%         cfg.subjectID            = subject(idx).id;
%         cfg.subjectPath          = fullfile(mPath, subject(idx).dir);
%         for sess = subject(idx).Restin
%             cfg.sessionID        = sess;
%             mnet_hcp_meeg_forward(cfg);
%         end
%     end
%     %% Run spDCM
%     cfg.TimeWin                  = [0 2000];
%     cfg.FreqWin                  = [1,100];
%     cfg.analysis                 = 'CSD';
%     cfg.model                    = 'ERP';
%     cfg.spatial                  = 'ECD';
%     % cfg.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27] [-28; -96; -6] [28;-96;-6]];
%     % cfg.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC', 'lV1','rV1'};
%     cfg.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27]];
%     cfg.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC'};
%     cfg.chanMode                 = 8;
% 
%     % cfg.A{1}                     = [0 0 0 0; 0 0 0 0; 1 1 0 0 ; 1 1 1 0];      % Forward
%     % cfg.A{2}                     = cfg.A{1}';                                  % Backward
%     % cfg.A{3}                     = [0 1 0 0 ;1 0 0 0; 0 0 0 0; 0 0 0 0];       % Modulatory
%     cfg.B{1}                     = zeros(4,4);
% 
%     cfg.A{1}                     = ones(4,4);
%     cfg.A{2}                     = ones(4,4);
%     cfg.A{3}                     = ones(4,4);
% 
%     cfg.effect                   = zeros(1,0); % Between trial effect (Design matrix)
%     cfg.effectName               = {};
% 
%     if ~isdir(fullfile(cfg.savePath,'DCMs'))
%         mkdir(fullfile(cfg.savePath,'DCMs'));
%     end
% 
%     for idx = 1:length(subject)
%         cfg.subjectID            = subject(idx).id;
%         cfg.subjectPath          = fullfile(mPath, subject(idx).dir);
%         for sess = subject(idx).Restin
%             cfg.sessionID        = sess;
%             cfg.trials           = [1]; % D.condlist{1}, LowAlpha
%             [D, DCM_Low]   = mnet_hcp_meeg_rest_spDCM_erpmodel(cfg);
%             save(fullfile(cfg.savePath,'DCMs',['DCM_' cfg.subjectID '_' num2str(cfg.sessionID) '_Low' upper(cfg.band(1)) cfg.band(2:end)]),'DCM_Low');
%             cfg.trials           = [2]; % D.condlist{2}, HighAlpha
%             [D, DCM_High]  = mnet_hcp_meeg_rest_spDCM_erpmodel(cfg);
%             save(fullfile(cfg.savePath,'DCMs',['DCM_' cfg.subjectID '_' num2str(cfg.sessionID) '_High' upper(cfg.band(1)) cfg.band(2:end)]),'DCM_High');
%             close all;
%         end
%     end
% end