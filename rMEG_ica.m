%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCP MEG ICA USING FAST ICA $$TEST VERSION$$                             %
%                                                                         %
% This code is for ica of HCP dataset                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.04.17 20:42 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Subject from HCP Dataset
addpath('/home/kuro/Codes/M_code/my_code');
addpath('/home/kuro/matlabwork/FastICA_25');
use_my_library('ft',1);
clear;
mPath='/media/kuro/DATA1/HCPMEG';
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

workingDir = '/media/kuro/DATA1/MEGResult';

%%
%% Set Path and Load Files
subjectPath = fullfile(mPath, subject(1).dir);
rMEGPath = fullfile(subjectPath,[subject(1).id '_MEG_Restin_preproc'],subject(1).id,'MEG','Restin','rmegpreproc');
rMEGRawData=[subject(1).id '_MEG_' num2str(subject(1).Restin(1)) '-Restin_rmegpreproc'];
% icaInfo = fullfile(subjectPath,[subject(1).id '_MEG_Restin_preproc'],subject(1).id,'MEG','Restin','icaclass',[subject(1).id,'_MEG_' num2str(subject(1).Restin(1)) '-Restin_icaclass_vs.mat']);
load(fullfile(rMEGPath,rMEGRawData)); % load HCP rMEG data
layoutfile = '4D248_helmet.mat'; % You need to set fieldtrip/template path.

anaPath = fullfile(subjectPath, [subject(1).id '_MEG_anatomy'],subject(1).id, 'MEG','anatomy');
headmodelfile = fullfile(anaPath,[subject(1).id '_MEG_anatomy_headmodel']);
sourcemodelfile = fullfile(anaPath,[subject(1).id '_MEG_anatomy_sourcemodel_2d']);
coordtransformfile = [subject(1).id '_MEG_anatomy_transform.txt'];

%%     Centering and Scaling for FastICA     %%
%% Demean and Scale every channels and every tials
data.demean = [];
data.mean = [];
data.scale = [];
data.demeanscaled = [];
data.trialsize = [];
data.concat = [];
for i = 1:length(data.trial)
    data.mean{i} = mean(data.trial{i},2);
    for j = 1:length(data.mean{i})
        data.demean{i}(j,:) = data.trial{i}(j,:) - data.mean{i}(j); % Demean every CHANNEL to have 0 mean
    end
    data.scale{i} = max(max(abs(data.demean{i}))); % Scale every TRIAL by its max absolute amplitude
    data.demeanscaled{i} = data.demean{i}/data.scale{i};
    data.trialsize{i} = size(data.trial{i} ,2);
    data.concat = [data.concat data.demeanscaled{i}];
end

[icasig A W]  = fastica(data.concat);


%%     SPECTRUM ANALYSIS    %%
%% Calculate the Powersepctrum
cfg = [];
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq =1;
cfg.keeptrials = 'no';
datapow = ft_freqanalysis(cfg, data);


%% plot powerspectrum
freqband = [2 4;5 7;8 12;15 30;30 60;60 90];
freqname = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma1', 'Gamme2'};

figure;
cfg = [];
cfg.layout = '4D248_helmet.mat';
for i =1:size(freqband, 1)
    cfg.xlim   = freqband(i,:);
    subplot(2,3,i);
    title([freqname{i} ' (' num2str(freqband(i,1)) '-' num2str(freqband(i,2)) 'Hz)']);
    ft_topoplotER(cfg, datapow);
    colorbar;
end

figure;
cfg = [];
cfg.channel = {'A129'}; % If multiple channel is set, then the function ft_singleplotER will calculate average power.
ft_singleplotER(cfg, datapow);

clear freqband freqname i;

%% Calculate the Powersepctrum of ic 1 signal
mixsig = A(:,1:150)*icasig(1:150,:);
reconSig = [];
for i = 1:147
    idxStart = ((i-1)*1018)+1;
    idxEnd = (i)*1018;
    reconSig{i}(:,:) = mixsig(:,idxStart:idxEnd);
end

for i = 1:length(reconSig)
    reconSig{i} = reconSig{i}*data.scale{i};
    for j = 1:length(data.mean{i})
        reconSig{i}(j,:) = reconSig{i}(j,:) + data.mean{i}(j);
    end
end
reconData = data;
for i = 1:length(reconSig)
    reconData.trial{i} = reconSig{i};
end

cfg = [];
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq =1;
cfg.keeptrials = 'no';
datapow = ft_freqanalysis(cfg, reconData);


%% plot powerspectrum
freqband = [2 4;5 7;8 12;15 30;30 60;60 90];
freqname = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma1', 'Gamme2'};

figure;
cfg = [];
cfg.layout = '4D248_helmet.mat';
for i =1:size(freqband, 1)
    cfg.xlim   = freqband(i,:);
    subplot(2,3,i);
    title([freqname{i} ' (' num2str(freqband(i,1)) '-' num2str(freqband(i,2)) 'Hz)']);
    ft_topoplotER(cfg, datapow);
    colorbar;
end

figure;
cfg = [];
cfg.channel = {'A129'}; % If multiple channel is set, then the function ft_singleplotER will calculate average power.
ft_singleplotER(cfg, datapow);

clear freqband freqname i;
