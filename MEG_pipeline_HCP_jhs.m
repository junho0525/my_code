%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCP MEG Source Reconstruction Pipeline Using SPM                        %
%                                                                         %
% This pipeline code is for HCP dataset                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.02.20 13:59 - By Junho Son                                     %
%     2018.03.14 06:44 - By Junho Son                                     %
%         - Add Fieldtrip Batch                                           %
%     2018.04.02 20:48 - By Junho Son                                     %
%         - Add Resting FT Batch(rMEG_pipeline_HCP_jhs_v4                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Subject from HCP Dataset
addpath('/home/kuro/Codes/M_code/my_code');
clear;
mPath='/projects1/HCPMEG';
megDirInfo = dir(mPath);
subject = [];
for i=1:length(megDirInfo)
    if ~isempty(strfind(megDirInfo(i).name, 'MEG')) % if subject folder is '[HCP data path]/[subject id]_MEG'
        subject(end+1).id = strtok(megDirInfo(i).name,'_');
        subject(end).dir = megDirInfo(i).name;
    elseif ~isempty(str2num(megDirInfo(i).name)) % if subject folder is '[HCP data path]/[subject id]'
        subject(end+1).id = megDirInfo(i).name;
        subject(end).dir = megDirInfo(i).name;
    end
end
clear megDirInfo;
for i =1:length(subject)
    subjDirInfo = dir(fullfile(mPath,subject(i).dir));
    subject(i).anatomy = [];
    subject(i).Restin = [];
    subject(i).Wrkmem = [];
    subject(i).StoryM = [];
    subject(i).Motort = [];
    for j =1:length(subjDirInfo)
       if ~isempty(strfind(subjDirInfo(j).name, 'anatomy'))
           subject(i).anatomy = 1; 
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Restin_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/Restin/rmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).Restin(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).Restin =unique(subject(i).Restin); 
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Wrkmem_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/Wrkmem/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).Wrkmem(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).Wrkmem =unique(subject(i).Wrkmem);
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'StoryM_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/StoryM/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).StoryM(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).StoryM =unique(subject(i).StoryM);
       end
       if ~isempty(strfind(subjDirInfo(j).name, 'Motort_preproc'))
           tempDirInfo = [];
           tempDirInfo = dir(fullfile(mPath, subject(i).dir,subjDirInfo(j).name,subject(i).id, 'MEG/Motort/tmegpreproc'));
           for k = 1:length(tempDirInfo)
               sessionNum =strsplit(strtok(tempDirInfo(k).name,'-'),'_MEG_');
               if length(sessionNum)>1&&~isempty(str2num(sessionNum{2}))
                   subject(i).Motort(end+1) = str2num(sessionNum{2}); 
               end
           end
           subject(i).Motort =unique(subject(i).Motort);
       end
    end
end
clear tempDirInfo subjDirInfo subDirInfo sessionNum i j k;

workingDir = '/projects1/MEGResult';

%% Example
% configuration for tMEG_pipeline_HCP_jhs_v3(cfg)
cfg = [];
% You must specify these five fields
cfg.taskType                = 'Wrkmem';
cfg.subjectID               = subject(8).id;
cfg.subjectPath             = fullfile(mPath,subject(8).dir);
cfg.savePath                = workingDir;
cfg.sessionID               = subject(8).Wrkmem(1);

% These fields are optional. If you do not specify, use defualt.
cfg.sourceLocalization      = 0; % yes : 1(default), no : 0 
cfg.showERFAll              = 0; % yes : 1, no : 0(defualt)
cfg.showERFHighSNR          = 0; % yes : 1(default), no : 0
cfg.showERFSNRThreshold     = 3; % We want to show ERF with high SNR. (defualt : 3)
cfg.showERFSNRMaxNum        = 10; % The numbers of ERFs that you want to show. (default : 20)
cfg.showForwardResult       = 0; % yes : 1, no : 0(defualt)
cfg.saveFig                 = 0; % yes : 1, no : 0(defualt)

% DCM parameters
cfg.DCM                     = 1;
cfg.DCMTimeWin              = [0 500];
cfg.DCManalysis             = 'ERP';
cfg.DCMROI = [ [-46;19;22] [-31;21;4] [-43;-76;35]  [-2;50.5;1.7]  [-28;-96;-6]];
cfg.DCMROIname = {'lDLPFC','lVLPRF' , 'lAG', 'LMPFC', 'lV1'};

cfg.DCMeffect = [1 0 ;0 1 ; 0 0];
cfg.DCMeffectName = {'0-back', '2-back'};
% cfg.DCMROI = [ [-46;19;22] [-31;21;4] [-43;-76;35] [-3;-54;31] [-2;50.5;1.7]  [-28;-96;-6]];
% cfg.DCMROIname = {'lDLPFC','lVLPRF' , 'lAG','lPcPr', 'LMPFC', 'lV1'};
% check ROIs
dip.n_seeds = 1;
dip.n_dip = size(cfg.DCMROI,2);
dip.Mtb = 1;
dip.j{1} = zeros(3*dip.n_dip, 1);
dip.loc{1} = cfg.DCMROI;
spm_eeg_inv_ecd_DrawDip('Init', dip);


DCM = tMEG_pipeline_HCP_jhs_v3(cfg);
%D{1}.con = 2;
%spm_eeg_invert_display(D{1});
%% Batch
cfg = [];
cfg.taskType = 'Wrkmem';
cfg.savePath = workingDir;
cfg.sourceLocalization      = 1; % yes : 1(default), no : 0 
cfg.showERFAll              = 0; % yes : 1, no : 0(defualt)
cfg.showERFHighSNR          = 1; % yes : 1(default), no : 0
cfg.showERFSNRThreshold     = 3; % We want to show ERF with high SNR. (defualt : 3)
cfg.showERFSNRMaxNum        = 10; % The numbers of ERFs that you want to show. (default : 20)
cfg.showForwardResult       = 0; % yes : 1, no : 0(defualt)
cfg.saveFig                 = 1; % yes : 1, no : 0(defualt)
cfg.DCM                     = 0;
cfg.DCManalysis             = 'CSD';

for i = 1:length(subject)
    cfg.subjectID           = subject(i).id;
    cfg.subjectPath         = fullfile(mPath, subject(i).dir);
    sess = [];
    for sess = subject(i).Wrkmem
        cfg.sessionID       = sess;
        D = tMEG_pipeline_HCP_jhs_v3(cfg);
        close all;
    end
end
%%
cfg = [];
cfg.taskType = 'Motort';
cfg.savePath = workingDir;
cfg.sourceLocalization      = 1; % yes : 1(default), no : 0 
cfg.showERFAll              = 0; % yes : 1, no : 0(defualt)
cfg.showERFHighSNR          = 1; % yes : 1(default), no : 0
cfg.showERFSNRThreshold     = 3; % We want to show ERF with high SNR. (defualt : 3)
cfg.showERFSNRMaxNum        = 10; % The numbers of ERFs that you want to show. (default : 20)
cfg.showForwardResult       = 0; % yes : 1, no : 0(defualt)
cfg.saveFig                 = 1; % yes : 1, no : 0(defualt)
cfg.DCM                     = 0;

for i = 1:length(subject)
    cfg.subjectID           = subject(i).id;
    cfg.subjectPath         = fullfile(mPath, subject(i).dir);
    for sess = subject(i).Motort
        cfg.sessionID       = sess;
        D = tMEG_pipeline_HCP_jhs_v3(cfg);
        close all;
    end
end

%% Fieldtrip ERP and Source Display
cfg = [];
cfg.taskType = 'Motort';
cfg.savePath = workingDir;
cfg.sourceLocalization      = 1; % yes : 1(default), no : 0 
cfg.showERFAll              = 1; % yes : 1, no : 0(defualt)
% cfg.showERFHighSNR          = 1; % yes : 1(default), no : 0
% cfg.showERFSNRThreshold     = 3; % We want to show ERF with high SNR. (defualt : 3)
% cfg.showERFSNRMaxNum        = 10; % The numbers of ERFs that you want to show. (default : 20)
cfg.showForwardResult       = 0; % yes : 1, no : 0(defualt)
cfg.saveFig                 = 1; % yes : 1, no : 0(defualt)
cfg.showMEGSignal           = 0;
cfg.showSensPower           = 1;

for i = 1:length(subject)
    cfg.subjectID           = subject(i).id;
    cfg.subjectPath         = fullfile(mPath, subject(i).dir);
    for sess = subject(i).Motort
        cfg.sessionID       = sess;
        D = tMEG_pipeline_HCP_jhs_v4(cfg);
        close all;
    end
end

%% Fieldtrip Resting Source Level FC Analysis
cd(workingDir);
restWorkingDir = fullfile(workingDir, 'Restin');
if ~isdir(restWorkingDir),mkdir(restWorkingDir);end

cfg = [];
cfg.taskType = 'Restin'; % Default
cfg.savePath = restWorkingDir;
cfg.showMEGSignal            = 0;
cfg.showForwardResult        = 0; % yes : 1, no : 0(default)
cfg.sourceLocalization       = 1;
cfg.showSensPower            = 0;
cfg.showAlphaPowRatio        = 0;
cfg.showConnectionDensityMap = 0;
cfg.showFullFC               = 0;
cfg.showSeedBasedFC          = 0; % yes : 1(default) no : 0
cfg.saveFig                  = 0;

for i = 1:length(subject)
    cfg.subjectID            = subject(i).id;
    cfg.subjectPath          = fullfile(mPath, subject(i).dir);
    for sess = subject(i).Restin
        cfg.sessionID        = sess;
        D = rMEG_pipeline_HCP_jhs_v4(cfg);
        close all;
    end
end
%% spDCM for Resting data
use_my_library('ft',0);
use_my_library('spm',1);
%
cd(workingDir);
restWorkingDir = fullfile(workingDir, 'Restin');
if ~isdir(restWorkingDir),mkdir(restWorkingDir);end

cfg                          = [];
cfg.taskType                 = 'Restin'; 
cfg.savePath                 = restWorkingDir;
cfg.showMEGSignal            = 0;        % default 0
cfg.showForwardResult        = 0;

cfg.TimeWin                  = [0 2000];
cfg.FreqWin                  = [1,100];
cfg.analysis                 = 'CSD';
cfg.model                    = 'ERP';
cfg.spatial                  = 'ECD';
cfg.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27]];
cfg.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC'};
cfg.trials                   = [1];
cfg.chanMode                 = 8;

cfg.A{1}                     = [0 0 0 0; 0 0 0 0; 1 1 0 0 ; 1 1 1 0];      % Forward
cfg.A{2}                     = cfg.A{1}';                                  % Backward
cfg.A{3}                     = [0 1 0 0 ;1 0 0 0; 0 0 0 0; 0 0 0 0];       % Modulatory
cfg.B{1}                     = zeros(4,4);

cfg.A{1}                     = randn(4,4)/12;
cfg.A{2}                     = randn(4,4)/12;
cfg.A{3}                     = randn(4,4)/12;

cfg.effect                   = zeros(1,0);
cfg.effectName               = {};

cfg.subjectID                = subject(1).id;
cfg.subjectPath              = fullfile(mPath, subject(1).dir);
cfg.sessionID                = 3;
[D, DCM]                     = rMEG_pipeline_HCP_jhs_v7(cfg);
%% spDCM for compare High&Low Alpha
use_my_library('ft',0);
use_my_library('spm',1);
%
cd(workingDir);
restWorkingDir = fullfile(workingDir, 'Restin');
if ~isdir(restWorkingDir),mkdir(restWorkingDir);end

cfg                          = [];
cfg.taskType                 = 'Restin'; 
cfg.savePath                 = restWorkingDir;
cfg.showMEGSignal            = 0;        % default 0
cfg.showForwardResult        = 0;

cfg.TimeWin                  = [0 2000];
cfg.FreqWin                  = [1,100];
cfg.analysis                 = 'CSD';
cfg.model                    = 'ERP';
cfg.spatial                  = 'ECD';
% cfg.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27] [-28; -96; -6] [28;-96;-6]];
% cfg.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC', 'lV1','rV1'};
cfg.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27]];
cfg.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC'};
cfg.chanMode                 = 8;

% cfg.A{1}                     = [0 0 0 0; 0 0 0 0; 1 1 0 0 ; 1 1 1 0];      % Forward
% cfg.A{2}                     = cfg.A{1}';                                  % Backward
% cfg.A{3}                     = [0 1 0 0 ;1 0 0 0; 0 0 0 0; 0 0 0 0];       % Modulatory
cfg.B{1}                     = zeros(4,4);

cfg.A{1}                     = ones(4,4);
cfg.A{2}                     = ones(4,4);
cfg.A{3}                     = ones(4,4);

cfg.effect                   = zeros(1,0);
cfg.effectName               = {}; % 'High Alpha'


for i = 1:length(subject)
    cfg.subjectID            = subject(i).id;
    cfg.subjectPath          = fullfile(mPath, subject(i).dir);
    for sess = subject(i).Restin
        cfg.sessionID        = sess;
        cfg.trials           = [1];
        [D, DCM_Low_Alpha]   = rMEG_pipeline_HCP_jhs_v7(cfg);
        save(fullfile(cfg.savePath,['DCM_' cfg.subjectID '_' num2str(cfg.sessionID) '_LowAlpha']),'DCM_Low_Alpha');
        cfg.trials           = [2];
        [D, DCM_High_Alpha]  = rMEG_pipeline_HCP_jhs_v7(cfg);
        save(fullfile(cfg.savePath,['DCM_' cfg.subjectID '_' num2str(cfg.sessionID) '_HighAlpha']),'DCM_High_Alpha');
        close all;
    end
end

%% DCM for further approximation
DCM2.name                    = DCM.name;
DCM2.xY.Dfile                = DCM.xY.Dfile;
DCM2.val                     = DCM.val;
DCM2.options                 = DCM.options;
DCM2.Lpos                    = DCM.Lpos;
DCM2.Sname                   = DCM.Sname;
DCM2.M                       = DCM.M;
DCM2.xU                      = DCM.xU;
DCM2.A                       = DCM.A;
DCM2.B                       = DCM.B;
DCM2.C                       = DCM.C;
DCM2.M.pE                    = DCM.Ep;
DCM                          = spm_dcm_csd_jhs(DCM2);

%% Check ROIs
dip.n_seeds = 1;
dip.n_dip = size(cfg.ROI,2);
dip.Mtb = 1;
dip.j{1} = zeros(3*dip.n_dip, 1);
dip.loc{1} = cfg.ROI;
spm_eeg_inv_ecd_DrawDip('Init', dip);

%% spDCM for compare source level High&Low band power
use_my_library('ft',0);
use_my_library('spm',1);
%
cd(workingDir);
restWorkingDir = fullfile(workingDir, 'Restin');
if ~isdir(restWorkingDir),mkdir(restWorkingDir);end

cfg                          = [];
cfg.taskType                 = 'Restin'; 
cfg.savePath                 = restWorkingDir;

cfg.showForwardResult        = 0;
cfg.band                     = 'alpha';
cfg.showPowerRatio           = 1;

cfg.TimeWin                  = [0 2000];
cfg.FreqWin                  = [1,100];
cfg.analysis                 = 'CSD';
cfg.model                    = 'ERP';
cfg.spatial                  = 'ECD';
% cfg.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27] [-28; -96; -6] [28;-96;-6]];
% cfg.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC', 'lV1','rV1'};
cfg.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27]];
cfg.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC'};
cfg.chanMode                 = 8;

% cfg.A{1}                     = [0 0 0 0; 0 0 0 0; 1 1 0 0 ; 1 1 1 0];      % Forward
% cfg.A{2}                     = cfg.A{1}';                                  % Backward
% cfg.A{3}                     = [0 1 0 0 ;1 0 0 0; 0 0 0 0; 0 0 0 0];       % Modulatory
cfg.B{1}                     = zeros(4,4);

cfg.A{1}                     = ones(4,4);
cfg.A{2}                     = ones(4,4);
cfg.A{3}                     = ones(4,4);

cfg.effect                   = zeros(1,0);
cfg.effectName               = {}; % 'High Alpha'


for i = 1:length(subject)
    cfg.subjectID            = subject(i).id;
    cfg.subjectPath          = fullfile(mPath, subject(i).dir);
    for sess = subject(i).Restin
        cfg.sessionID        = sess;
        cfg.trials           = [1];
        [D, DCM_Low_Alpha]   = rMEG_pipeline_HCP_jhs_v8(cfg);
        save(fullfile(cfg.savePath,['DCM_' cfg.subjectID '_' num2str(cfg.sessionID) '_LowAlpha']),'DCM_Low_Alpha');
        cfg.trials           = [2];
        [D, DCM_High_Alpha]  = rMEG_pipeline_HCP_jhs_v8(cfg);
        save(fullfile(cfg.savePath,['DCM_' cfg.subjectID '_' num2str(cfg.sessionID) '_HighAlpha']),'DCM_High_Alpha');
        close all;
    end
end

%% Sensor/Source level analysis for compare source level High&Low Band Power
use_my_library('ft',0);
use_my_library('spm',1);
%
cd(workingDir);
restWorkingDir = fullfile(workingDir, 'Restin');
if ~isdir(restWorkingDir),mkdir(restWorkingDir);end

cfg                          = [];
cfg.taskType                 = 'Restin'; 
cfg.savePath                 = restWorkingDir;

cfg.showForwardResult        = 0;
cfg.band                     = 'alpha';
cfg.showPowerRatio           = 1;

cfg.TimeWin                  = [0 2000];
cfg.FreqWin                  = [1,100];
cfg.analysis                 = 'CSD';
cfg.model                    = 'ERP';
cfg.spatial                  = 'ECD';
% cfg.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27] [-28; -96; -6] [28;-96;-6]];
% cfg.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC', 'lV1','rV1'};
cfg.ROI                      = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27]];
cfg.ROIName                  = {'lLP', 'rLP', 'Prec', 'mPFC'};
cfg.chanMode                 = 8;

% cfg.A{1}                     = [0 0 0 0; 0 0 0 0; 1 1 0 0 ; 1 1 1 0];      % Forward
% cfg.A{2}                     = cfg.A{1}';                                  % Backward
% cfg.A{3}                     = [0 1 0 0 ;1 0 0 0; 0 0 0 0; 0 0 0 0];       % Modulatory
cfg.B{1}                     = zeros(4,4);

cfg.A{1}                     = ones(4,4);
cfg.A{2}                     = ones(4,4);
cfg.A{3}                     = ones(4,4);

cfg.effect                   = zeros(1,0);
cfg.effectName               = {}; % 'High Alpha'


for i = 1:length(subject)
    cfg.subjectID            = subject(i).id;
    cfg.subjectPath          = fullfile(mPath, subject(i).dir);
    for sess = subject(i).Restin
        cfg.sessionID        = sess;
        D                    = power_ratio_map_visual_chan(cfg);
        % save(fullfile(cfg.savePath,['DCM_' cfg.subjectID '_' num2str(cfg.sessionID) '_LowAlpha']),'DCM_Low_Alpha');
    end
end
