function D = tMEG_pipeline_HCP_jhs_v4(config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCP Task MEG Source Reconstruction Pipeline Using FieldTrip             %
%                                                                         %
% This pipeline code is for task HCP dataset                              %
%                                                                         %
% D = tMEG_pipeline_HCP_jhs_v3(config)                                    %
%    Output      D     - spm meeg object                                  %
%    Input                                                                %
%       config.subjectID     - HCP dataset subject id                     %
%       config.tasktype      - 'Wrkmem', 'Motort', 'StoryM'               %
%       config.subjectPath   - HCP Single Subject Data Path               %
%       config.savePath      - Save path for output and intermediate files%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.03.14 09:55 - By Junho Son                                     %
%         - Make this code                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input
% These five fields must be specified
if ~isfield(config, 'subjectID'),error('Invlid Subject ID!! There is no subjectID');end
if ~isfield(config, 'taskType'),error('Invalid Task Type!! There is no taskType');end
if (~isfield(config, 'subjectPath')||~isdir(config.subjectPath));error('Invalid Subject Path!! There is no such folder');end
if (~isfield(config, 'savePath')||~isdir(config.savePath)),error('Invalid Save Path!! There is no such folder');end
if (~isfield(config, 'sessionID')),error('Invalid Session ID!! There is no session ID field!!');end

% These fields is optional. If these fields are not specified, use default.
if ~isfield(config, 'sourceLocalization'), config.sourceLocalization = 1; end % defulat : 1 
if ~isfield(config, 'showERFAll'), config.showERFAll = 0; end % defualt : 0
if ~isfield(config, 'showForwardResult'), config.showForwardResult = 0; end % defualt : 0
%if ~isfield(config, 'showERFHighSNR'), config.showERFHighSNR = 1; end
%if ~isfield(config, 'showERFSNRThreshold'), config.showERFSNRThreshold = 3; end
if ~isfield(config, 'showMEGSignal'), config.showMEGSignal = 0; end %%%%%
if ~isfield(config, 'saveFig'), config.saveFig = 0; end %%%%%
if ~isfield(config, 'showSensPower'), config.showSensPower = 0; end %%%%%
%% Set Path
switch config.taskType
    case 'Wrkmem'
        tMEGPath = fullfile(config.subjectPath,[config.subjectID '_MEG_Wrkmem_preproc'],config.subjectID,'MEG','Wrkmem','tmegpreproc');
        tMEGRawData=[config.subjectID '_MEG_' num2str(config.sessionID) '-Wrkmem_tmegpreproc_TIM'];
        icaInfo = fullfile(config.subjectPath,[config.subjectID '_MEG_Wrkmem_preproc'],config.subjectID,'MEG','Wrkmem','icaclass',[config.subjectID,'_MEG_' num2str(config.sessionID) '-Wrkmem_icaclass_vs.mat']);
    case 'Motort'
        tMEGPath = fullfile(config.subjectPath,[config.subjectID '_MEG_Motort_preproc'],config.subjectID,'MEG','Motort','tmegpreproc');
        tMEGRawData=[config.subjectID '_MEG_' num2str(config.sessionID) '-Motort_tmegpreproc_TFLA'];
        icaInfo = fullfile(config.subjectPath,[config.subjectID '_MEG_Motort_preproc'],config.subjectID,'MEG','Motort','icaclass',[config.subjectID,'_MEG_' num2str(config.sessionID) '-Motort_icaclass_vs.mat']);
    case 'StoryM'
        tMEGPath = fullfile(config.subjectPath,[config.subjectID '_MEG_StoryM_preproc'],config.subjectID,'MEG','StoryM','tmegpreproc');
        tMEGRawData=[config.subjectID '_MEG_' num2str(config.sessionID) '-StoryM_tmegpreproc_TIM'];
        icaInfo = fullfile(config.subjectPath,[config.subjectID '_MEG_StoryM_preproc'],config.subjectID,'MEG','StoryM','icaclass',[config.subjectID,'_MEG_' num2str(config.sessionID) '-StoryM_icaclass_vs.mat']);
    otherwise
        Error();
end
anaPath = fullfile(config.subjectPath, [config.subjectID '_MEG_anatomy'],config.subjectID, 'MEG','anatomy');
headmodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_headmodel']);
sourcemodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_sourcemodel_2d']);
layoutfile = '4D248_helmet.mat'; % You need to set fieldtrip/template path.
% Set Save Path
if config.saveFig
    mkdir(fullfile(config.savePath, 'figure'));
    mkdir(fullfile(config.savePath, 'figure',config.subjectID));
    mkdir(fullfile(config.savePath, 'figure',config.subjectID, config.taskType));
    savePath = fullfile(config.savePath, 'figure',config.subjectID, config.taskType);
end

%%     PREPROCESSING     %%
%% Load MEG data
load(fullfile(tMEGPath,tMEGRawData)); % load HCP rMEG data
cfg                        = [];
cfg.trials                 = 'all';
switch config.taskType
    case 'Wrkmem'
        cfg.toilim         = [-0.5 1.5];
    case 'Motort'
        cfg.toilim         = [-0.3 1.2];
    otherwise
        cfg.toilim         = [-0.5 2];
end
data                       = ft_redefinetrial(cfg, data) % redefine timewindow of trial 
%% Baseline correction
cfg.demean                 = 'yes';
cfg.detrend                = 'yes';
cfg.baselinewindow         = [-0.2 0];
cfg.lpfilter               = 'yes';
cfg.lpfreq                 = 30;
cfg.hpfilter               = 'yes';
cfg.hpfreq                 = 1;
data                       = ft_preprocessing(cfg, data);
%% Browse data.
if config.showMEGSignal
    cfg                    = []; 
    cfg.layout             = layoutfile;
    cfg.continuous         = 'no';
    cfg.viewmode           = 'vertical';
    ft_databrowser(cfg, data);
end
%%     SPECTRUM ANALYSIS    %%
if config.showSensPower
    %% Calculate the Powersepctrum
    cfg = [];
    cfg.output = 'pow';
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.tapsmofrq =1;
    cfg.keeptrials = 'no';
    datapow = ft_freqanalysis(cfg, data);
    %% plot powerspectrum
    figure('units','normalized','outerposition', [0,0,1,1]);
    cfg = [];
    cfg.layout = '4D248_helmet.mat';
    cfg.xlim   = [1 30]; % frequency band
    title('Powerspectrum of 1 - 30 Hz');
    ft_topoplotER(cfg, datapow);
    colorbar;
    if config.saveFig
        mkdir(fullfile(savePath, 'SensPower'));
        saveas(gcf ,fullfile(savePath, 'SensPower',['FT_' config.subjectID '_' num2str(config.sessionID) '_sensPow.png']));
    end
end
%%     ERF ANALYSIS     %%
%% Make ERP for each conditions
% data.trialinfo - Trial information for TIM 
% data.trialinfo(:,4) - 1: face, 2: tools 0: fixation
% data.trialinfo(:,5) - 1: 0-back 2: 2-back NaN: fixation
% We define conditions only 0 back and 2 back for now.
cfg = [];
cfg.channel = 'MEG';
switch config.taskType
    case 'Wrkmem'
        cfg.trials = find(data.trialinfo(:,5)==1)'; % 0-back
        ERFdata(1) = ft_timelockanalysis(cfg, data);
        cfg.trials = find(data.trialinfo(:,5)==2)'; % 2-back
        ERFdata(2) = ft_timelockanalysis(cfg, data);
        cfg.trials = setdiff(setdiff(1:length(data.trial), find(data.trialinfo(:,5)==1)'),find(data.trialinfo(:,5)==2)'); % Fixation
        ERFdata(3) = ft_timelockanalysis(cfg, data);
    case 'Motort'
        cfg.trials = find(data.trialinfo(:,2)==1)'; % Left hand
        ERFdata(1) = ft_timelockanalysis(cfg, data);
        cfg.trials = find(data.trialinfo(:,2)==2)'; % Left foot
        ERFdata(2) = ft_timelockanalysis(cfg, data);
        cfg.trials = find(data.trialinfo(:,2)==4)'; % Right hand
        ERFdata(3) = ft_timelockanalysis(cfg, data);
        cfg.trials = find(data.trialinfo(:,2)==5)'; % Right foot
        ERFdata(4) = ft_timelockanalysis(cfg, data);
        cfg.trials = find(data.trialinfo(:,2)==6)'; % Fixation
        ERFdata(5) = ft_timelockanalysis(cfg, data);
        ERFdata(5).time = ERFdata(1).time;
        ERFdata(5).avg = [zeros(size(ERFdata(5).avg,1),size(ERFdata(1).avg,2)-size(ERFdata(5).avg,2)) ERFdata(5).avg]; % zero expansion with undefined fixation time point
    otherwise
        error();
end
%% Plot data on 2D Channel Layout
if config.showERFAll
    figure('units','normalized','outerposition', [0,0,1,1]);
    cfg                        = [];
    cfg.showlabels             = 'yes';
    cfg.showoutline            = 'yes';
    cfg.layout                 = layoutfile;
    switch config.taskType
        case 'Wrkmem'
            ft_multiplotER(cfg, ERFdata(1), ERFdata(2), ERFdata(3));
        case 'Motort'
            ft_multiplotER(cfg, ERFdata(1), ERFdata(2), ERFdata(3), ERFdata(4), ERFdata(5));
    end
    if config.saveFig
        mkdir(fullfile(savePath, 'ERF'));
        saveas(gcf, fullfile(savePath, 'ERF', ['FT_' config.subjectID '_' num2str(config.sessionID) '_ERF.png'] ))
    end
end
