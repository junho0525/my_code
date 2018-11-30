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
%     2018.05.14 12:12 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Subject from HCP Dataset
addpath('E:\my_code');
clear;
clear;
use_my_library init;
use_my_library('spm',1); % We are going to use spm.


mPath='E:\HCPMEG'; % This is the directory of root of HCP dataset.
workingDir = 'E:\MEGResult'; % This is the directory that all result files will be placed.
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

cd(workingDir);
restWorkingDir               = fullfile(workingDir, 'Restin');
if ~isdir(restWorkingDir),mkdir(restWorkingDir);end
dataset = cell(length(subject), 3, 2);
cfg.savePath                 = restWorkingDir;
%% Load and Do Preprocess
cfg.freqband                 = [0.1, 40];
cfg.downsample               = 200;

for i = 1:length(subject)
    cfg.subjectID            = subject(i).id;
    cfg.subjectPath          = fullfile(mPath, subject(i).dir);
    for sess = 1:length(subject(i).Restin)
        cfg.sessionID        = subject(i).Restin(sess);
        [dataset{i,sess,1},dataset{i,sess,2}] = mnet_hcp_meeg_preprocess_spm(cfg);
    end
end
% You should check by these two fuctions
% D = spm_eeg_load('filename');
% spm_eeg_review(D);
%% Set Forward Model
for i = 1:length(subject)
    cfg.subjectID            = subject(i).id;
    cfg.subjectPath          = fullfile(mPath, subject(i).dir);
    for sess = length(subject(i).Restin)
        cfg.sessionID        = subject(i).Restin(sess);
        [dataset{i,sess,1},dataset{i,sess,2}] = mnet_hcp_meeg_forward(cfg,dataset{i,sess,1},dataset{i,sess,2});
    end
end
%% Show Band-Contrast Source Activity (optional)
try
    cfg                      = rmfield(cfg, 'freqband');
    cfg                      = rmfield(cfg, 'downsample');
end
cfg.band                     = 'theta';
source_new = [];
source_mean = [];
for i = 1:length(subject)
    cfg.subjectID            = subject(i).id;
    cfg.subjectPath          = fullfile(mPath, subject(i).dir);
    for sess = 1:length(subject(i).Restin)
        cfg.sessionID        = subject(i).Restin(sess);
        source_new = [];
        source_new = mnet_plot_band_contrast(cfg,dataset{i,sess,2});
        if isempty(source_mean)
            source_mean = source_new;
        else
            source_mean.pow = source_new.pow + source_mean.pow;
        end
    end
end
bnd.pnt = sourcemodel2d.pos;
bnd.tri = sourcemodel2d.tri;
mask = ~~(source_sim);
alpha = 0.9;
mask = mask + (~mask)*alpha;
ft_plot_mesh(bnd);
ft_plot_mesh(bnd, 'vertexcolor', source_mean.pow, 'facealpha',source_mean.mask,'colormap','jet');
lighting gouraud
camlight

title('Simulated Source Location');
%% Redefine Trials According to Band Power
% set band : 'delta' (2-4Hz), 'theta' (5-7Hz), 'alpha' (8-12Hz), 
%    'beta' (15-30Hz), 'gamma1' (30-60Hz), 'gamma2(60-90Hz)
% Trials defined by median power of median frequency of each band. 
%    ex. for alpha, calculate 10Hz(median frequency of alpha) power for 
%        each trials and separate trials two group with larger than median
%        trial power and smaller than median trial power.
cfg.band                     = 'alpha';
cfg.showPowerRatio           = 0;
for i = 1:length(subject)
    cfg.subjectID            = subject(i).id;
    cfg.subjectPath          = fullfile(mPath, subject(i).dir);
    for sess = length(subject(i).Restin)
        cfg.sessionID        = subject(i).Restin(sess);
        dataset{i,sess,1}=mnet_hcp_meeg_definetrial_bandpower(cfg,dataset{i,sess,1},dataset{i,sess,2});
    end
end
%% Run spDCM
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

cfg.effect                   = zeros(1,0); % Between trial effect (Design matrix)
cfg.effectName               = {};

if ~isdir(fullfile(cfg.savePath,'DCMs'))
    mkdir(fullfile(cfg.savePath,'DCMs'));
end

for i = 1:length(subject)
    cfg.subjectID            = subject(i).id;
    cfg.subjectPath          = fullfile(mPath, subject(i).dir);
    for sess = subject(i).Restin
        cfg.sessionID        = sess;
        cfg.trials           = [1]; % D.condlist{1}, LowAlpha
        [D, DCM_Low]   = mnet_hcp_meeg_spDCM_erp_rest(cfg);
        save(fullfile(cfg.savePath,'DCMs',['DCM_' cfg.subjectID '_' num2str(cfg.sessionID) '_Low' upper(cfg.band(1)) cfg.band(2:end)]),'DCM_Low');
        cfg.trials           = [2]; % D.condlist{2}, HighAlpha
        [D, DCM_High]  = mnet_hcp_meeg_spDCM_erp_rest(cfg);
        save(fullfile(cfg.savePath,'DCMs',['DCM_' cfg.subjectID '_' num2str(cfg.sessionID) '_High' upper(cfg.band(1)) cfg.band(2:end)]),'DCM_High');
        close all;
    end
end
%% PEB
% load DCMs
dirInfo = dir(fullfile(restWorkingDir,'DCMs','Theta'));
fileRemoveIdx = zeros(0,1);
for i =1:length(dirInfo)
    fileRemoveIdx(end+1) = isempty(strfind(dirInfo(i).name, 'DCM'))|| isempty(strfind(dirInfo(i).name, [upper(cfg.band(1)) cfg.band(2:end)]));
end
dirInfo(find(fileRemoveIdx)) = []; % get filenames contain string 'DCM'.
clear fileRemoveIdx i

Descriptions = [];                                                         
DCMs = cell(0,2);                                                          

i = 1;
while i < length(dirInfo)
    fileNameSplit1 = strsplit(dirInfo(i).name,'_');
    fileNameSplit2 = strsplit(dirInfo(i+1).name,'_');
    if strcmp(fileNameSplit1{2},fileNameSplit2{2}) && strcmp(fileNameSplit1{3},fileNameSplit2{3})
        Descriptions(end+1).subjectId = fileNameSplit1{2};
        Descriptions(end).sessionId = fileNameSplit1{3};
        load(fullfile(fullfile(restWorkingDir,'DCMs','Theta'),dirInfo(i).name));
        load(fullfile(fullfile(restWorkingDir,'DCMs','Theta'),dirInfo(i+1).name));
        DCMs(end+1,:) = {DCM_High,DCM_Low};
    else
        i=i+1;
        continue;
    end
    i = i+2;
end
% Sort DCMs By sessionId
Descriptions_new= [];   
DCMs_old = DCMs;
DCMs = cell(0,6);
i = 1;
while i < length(Descriptions)
    if~(strcmp(Descriptions(i).subjectId, Descriptions(i+1).subjectId)&&...
            strcmp(Descriptions(i+1).subjectId,Descriptions(i+2).subjectId))
        i=i+1;
        continue;
    end
    Descriptions_new(end+1).subjectId =Descriptions(i).subjectId;
    DCMs(end+1,:) = {DCMs_old{i,1},DCMs_old{i+1,1},DCMs_old{i+2,1},DCMs_old{i,2},DCMs_old{i+1,2},DCMs_old{i+2,2}};
    i = i+3;
end
Descriptions = Descriptions_new;
clear DCMs_old DCM_High DCM_Low Descriptions_new i dirInfo fileNameSplit1 fileNameSplit2

PEBi = {}; M1.X = ones(6,2); M1.X(4:6,2) = -1;
for s = 1:size(DCMs, 1)
    dcm = {};
    for j=1:6, dcm{j,1}=DCMs{s,j};end
    PEBi{s,1} = spm_dcm_peb(dcm,M1);
end
use_my_library('mnet',1);
PEBa = spm_dcm_peb_peb(PEBi);
mnet_dcm_peb_display(PEBa);
use_my_library('mnet',0);
%%
freqbands = {'theta', 'delta', 'beta','gamma1','gamma2','alpha'};
for j = 4:6
    %% Redefine Trials According to Band Power
    try
        cfg                      = rmfield(cfg, 'freqband');
        cfg                      = rmfield(cfg, 'downsample');
    end
    % set band : 'delta' (2-4Hz), 'theta' (5-7Hz), 'alpha' (8-12Hz), 
    %    'beta' (15-30Hz), 'gamma1' (30-60Hz), 'gamma2(60-90Hz)
    % Trials defined by median power of median frequency of each band. 
    %    ex. for alpha, calculate 10Hz(median frequency of alpha) power for 
    %        each trials and separate trials two group with larger than median
    %        trial power and smaller than median trial power.
    cfg.band                     = freqbands{j};
    cfg.showPowerRatio           = 0;
    for i = 1:length(subject)
        cfg.subjectID            = subject(i).id;
        cfg.subjectPath          = fullfile(mPath, subject(i).dir);
        for sess = subject(i).Restin
            cfg.sessionID        = sess;
            mnet_hcp_meeg_definetrial_bandpower(cfg);
        end
    end
    %% Set Forward Model
    for i = 1:length(subject)
        cfg.subjectID            = subject(i).id;
        cfg.subjectPath          = fullfile(mPath, subject(i).dir);
        for sess = subject(i).Restin
            cfg.sessionID        = sess;
            mnet_hcp_meeg_forward(cfg);
        end
    end
    %% Run spDCM
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

    cfg.effect                   = zeros(1,0); % Between trial effect (Design matrix)
    cfg.effectName               = {};

    if ~isdir(fullfile(cfg.savePath,'DCMs'))
        mkdir(fullfile(cfg.savePath,'DCMs'));
    end

    for i = 1:length(subject)
        cfg.subjectID            = subject(i).id;
        cfg.subjectPath          = fullfile(mPath, subject(i).dir);
        for sess = subject(i).Restin
            cfg.sessionID        = sess;
            cfg.trials           = [1]; % D.condlist{1}, LowAlpha
            [D, DCM_Low]   = mnet_hcp_meeg_spDCM_erp_rest(cfg);
            save(fullfile(cfg.savePath,'DCMs',['DCM_' cfg.subjectID '_' num2str(cfg.sessionID) '_Low' upper(cfg.band(1)) cfg.band(2:end)]),'DCM_Low');
            cfg.trials           = [2]; % D.condlist{2}, HighAlpha
            [D, DCM_High]  = mnet_hcp_meeg_spDCM_erp_rest(cfg);
            save(fullfile(cfg.savePath,'DCMs',['DCM_' cfg.subjectID '_' num2str(cfg.sessionID) '_High' upper(cfg.band(1)) cfg.band(2:end)]),'DCM_High');
            close all;
        end
    end
end