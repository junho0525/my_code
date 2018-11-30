%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          This Is M/EEG task Source Localization Pipeline Code           %
%                           By Using Fieldtrip                            %
%                                                                         %
%   We need five inputs                                                   %
%   You need to specify these parameters.( All of them are string )       %
%       filetype        - File type               ('EEG', 'MEG' or 'HCP') %
%       rawfile         - M/EEG raw data filename (eg. Subject2Rest.ds)   %
%       subjectid       - Subject ID for HCP data (eg. '100307')          %
%       sourcemodelfile - Source model filename   (eg. sourcemodel4k.mat) %
%       headmodelfile   - Head model filename     (eg. headmodel.mat)     %
%       layoutfile      - Channel layout file     (eg. 4D248_helmet.mat)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
% -- 2018.02.13 xx:xx -- By Junho Son                                     %
% -- 2018.02.12 xx:xx -- By Junho Son                                     %
% -- 2018.02.28 15:16 -- By Junho Son                                     %
%                Redefine inputs and options as cfg struct                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    SET DATA PATH    %%
%% Check Input Parameters
%  If you want to use example datapath, please excute this section twice.
if ~exist('filetype')
   warning('You should specify file type.');
   filetype = 'HCP';
end
if exist('rawfile')
    israwfile=1;
else
    israwfile =0;
end
if exist('subjectid')
    issubjectid = 1;
else
    issubjectid = 0;
end

% distract path from rawfile
if ~israwfile&&~issubjectid
    warning('There is no rawfile. This code use default dataset (HCP-subject 100307 resting data)');
%     rawfile = '100307_MEG_7-Wrkmem_tmegpreproc_TIM.mat';
    rawfile = '/media/kuro/DATA1/HCPMEG/100307_MEG/100307_MEG_Wrkmem_preproc/100307/MEG/Wrkmem/tmegpreproc/100307_MEG_7-Wrkmem_tmegpreproc_TIM.mat';
%     megrootpath = '';
%     anapath = '';
%     tmegpath = '';
%     headmodelfile = '100307_MEG_anatomy_headmodel.mat';
%     sourcemodelfile = '100307_MEG_anatomy_sourcemodel_2d.mat';
end
switch filetype
    case 'HCP'
        if israwfile
            [tmegpath rawfile]    = fileparts(rawfile);
            subjectid = strtok(rawfile, '_');
            strind = strfind(tmegpath, subjectid);
            
            megrootpath = tmegpath(1:(strind(1)-2)); % '/media/kuro/DATA1/HCPMEG'
            tmegpath = tmegpath(strind(1):end); % '100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc'
            anapath = fullfile([subjectid, '_MEG'],[subjectid, '_MEG_anatomy'],subjectid, 'MEG','anatomy'); %100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy
            headmodelfile = [subjectid '_MEG_anatomy_headmodel.mat'];
            sourcemodelfile = [subjectid '_MEG_anatomy_sourcemodel_2d.mat'];
            
            clear strind;
        elseif issubjectid
            megrootpath = '/media/kuro/DATA1/HCPMEG'; % Default HCP MEG root dicrectory
            tmegpath = fullfile([subjectid, '_MEG'],[subjectid, '_MEG_Wrkmem_preproc'],subjectid,'MEG','Wrkmem','tmegpreproc');
            anapath = fullfile([subjectid, '_MEG'],[subjectid, '_MEG_anatomy'],subjectid, 'MEG','anatomy'); %100307_MEG/100307_MEG_anatomy/100307/MEG/anatomy
            headmodelfile = [subjectid '_MEG_anatomy_headmodel.mat'];
            sourcemodelfile = [subjectid '_MEG_anatomy_sourcemodel_2d.mat'];
            rawfile = [subfectid '_MEG_7-Wrkmem_tmegpreproc_TIM.mat']; 
        end
    case 'MEG' % This is for future work
        error('pipeline for general MEG data is not ready yet');
    case 'EEG' % This is for future work
        error('pipeline for general EEG data is not ready yet');
    otherwise
        error('unexpected filetype. filetype should be one of ''EEG'', ''MEG'', or ''HCP''.');
end
clear israwfile issubjectid filetype;

%%     PREPROCESSING     %%
%% Load MEG data

load(fullfile(megrootpath, tmegpath, rawfile)); % load HCP rMEG data
cfg                        = [];
cfg.trials                 = 'all';
cfg.toilim                 = [-0.5 2];
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
cfg                        = []; 
layoutfile                 = '4D248_helmet.mat'; % You need to set fieldtrip/template path.
[layoutpath layoutfile]    = fileparts(layoutfile); 
cfg.layout                 = layoutfile;
cfg.continuous             = 'no';
cfg.viewmode               = 'vertical';
ft_databrowser(cfg, data);

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
figure;
cfg = [];
cfg.layout = '4D248_helmet.mat';
cfg.xlim   = [1 30]; % frequency band
title('Powerspectrum of 1 - 30 Hz');
ft_topoplotER(cfg, datapow);
colorbar;

figure;
cfg = [];
cfg.channel = {'A129'}; % If multiple channel is set, then the function ft_singleplotER will calculate average power.
ft_singleplotER(cfg, datapow);

%%     ERF ANALYSIS     %%
%% Make ERP for each conditions
% data.trialinfo - Trial information for TIM 
% data.trialinfo(:,4) - 1: face, 2: tools 0: fixation
% data.trialinfo(:,5) - 1: 0-back 2: 2-back NaN: fixation
% We define conditions only 0 back and 2 back for now.
cfg = [];
cfg.channel = 'MEG';
cfg.trials = find(data.trialinfo(:,5)==1)'; % 0-back
zerobackdata = ft_timelockanalysis(cfg, data);

cfg.trials = find(data.trialinfo(:,5)==2)'; % 2-back
twobackdata = ft_timelockanalysis(cfg, data);

cfg.trials = setdiff(setdiff(1:length(data.trial), find(data.trialinfo(:,5)==1)'),find(data.trialinfo(:,5)==2)');
fixdata = ft_timelockanalysis(cfg, data);

%% Plot erp data on 2D channel layout
cfg                        = [];
cfg.showlabels             = 'yes';
cfg.showoutline            = 'yes';
cfg.layout                 = layoutfile;
cfg.baseline               = [-0.2,0];
cfg.linewidth              = 1;
cfg.xlim                   = [-0.2 1.5];
cfg.ylim                   = [-3e-13 3e-13];
figure('units', 'normalized', 'outerposition',[0,0,1,1]); % full screen window
ft_multiplotER(cfg, zerobackdata, twobackdata, fixdata);
%% Plot erp data as topoplot
cfg                        = [];
cfg.layout                 = layoutfile;
cfg.xlim                   = [-0.2:0.1:1.0]; % Define 12 time intervals
cfg.zlim                   = [-1e-13 1e-13]; % Set the color limits
figure;
ft_topoplotER(cfg, zerobackdata);
suptitle('0-back');
figure;
ft_topoplotER(cfg, twobackdata);
suptitle('2-back');
figure
ft_topoplotER(cfg, fixdata);
suptitle('Fixation');
%% Plot Channels According to SNR ratio
for i = 1:size(data.label,1)
    r(i) = snr(zerobackdata.avg(i,:), fixdata.avg(i,:));
end
ind                        = find(r>3); % ind : channels that snr ratio is greater than 3.
if size(ind)>20
    subplotnum             = 20;
else
    subplotnum             = length(ind);
end
figure;
for j = 1:ceil(subplotnum/5)
    for i = 1:5
        idx = (j-1)*4+i;
        if idx>length(ind)
            break
        end
        subplot(ceil(subplotnum/5),5,idx);
        
        hold on;
        plot(fixdata.time, fixdata.avg(ind(idx),:), 'c');
        plot(zerobackdata.time, zerobackdata.avg(ind(idx),:),'b');
        plot(twobackdata.time, twobackdata.avg(ind(idx),:), 'r');
        title(data.label(ind(idx)));
        hold off;
        xlim([-0.5 2]);
        ylim([-3e-13 3e-13]);
    end
end
clear i j ind idx r subplotnum;

%% Calculate the planar gradient
cfg                        = [];
cfg.feedback               = 'no'; % Show Neighbours of each channels
cfg.method                 = 'template';
cfg.neighbours             = ft_prepare_neighbours(cfg, zerobackdata);
cfg.planarmethod           = 'sincos';
zerobackdataplanar         = ft_megplanar(cfg, zerobackdata);

cfg.neighbours             = ft_prepare_neighbours(cfg, twobackdata);
twobackdataplanar          = ft_megplanar(cfg, twobackdata);

cfg.neighbours             = ft_prepare_neighbours(cfg, fixdata);
fixdataplanar          = ft_megplanar(cfg, fixdata);

cfg                        = [];
zerobackdataplanarComb     = ft_combineplanar(cfg, zerobackdataplanar);
twobackdataplanarComb     = ft_combineplanar(cfg, twobackdataplanar);
fixdataplanarComb     = ft_combineplanar(cfg, fixdataplanar);

%% Plot the results (planar gradients)
cfg                        = [];
cfg.layout                 = layoutfile;
cfg.xlim                   = [-0.2:0.1:1.0]; % Define 12 time intervals
cfg.zlim                   = [-1.2e-13 1.2e-13]; % Set the color limits
cfg.colorbar               = 'no';
figure;
colormap('jet');
ft_topoplotER(cfg, zerobackdata);
suptitle('0-Back Axial Gradient');
figure;
colormap('jet');
ft_topoplotER(cfg, twobackdata);
suptitle('2-Back Axial Gradient');
figure;
colormap('jet');
ft_topoplotER(cfg, fixdata);
suptitle('Fixation Axial Gradient');
cfg.zlim                   = [-1.25e-11 1.25e-11]; % Set the color limits
figure;
colormap('jet');
ft_topoplotER(cfg, zerobackdataplanarComb);
suptitle('0-Back Planar Gradient');
figure;
colormap('jet');
ft_topoplotER(cfg, twobackdataplanarComb);
suptitle('2-Back Planar Gradient');
figure;
colormap('jet');
ft_topoplotER(cfg, fixdataplanarComb);
suptitle('Fixation Planar Gradient');
%%    SOURCE RECONSTRUCTION OF ERP    %%
%% Forward solution(Calculate leadfield)
headmodel = load(fullfile(megrootpath, anapath, headmodelfile));
headmodel = headmodel.headmodel;
sourcemodel = load(fullfile(megrootpath, anapath ,sourcemodelfile));
sourcemodel = sourcemodel.sourcemodel2d;

cfg                        = [];
cfg.grad                   = zerobackdata.grad; % sensor position
cfg.channel                = {'MEG'};
cfg.grid                   = sourcemodel; % sourcemodel
cfg.headmodel              = headmodel; % volume conduction model
leadfield                  = ft_prepare_leadfield(cfg); % We can calculate leadfield without data itself.

clear headmodelfile sourcemodelfile;
%% Source Estimation
cfg                        = [];
cfg.method                 = 'mne'; % Minimum Norm Estimation
cfg.grid                   = leadfield;
cfg.headmodel              = headmodel;
cfg.mne.prewhiten          = 'yes';
cfg.mne.lambda             = 3;
cfg.mne.scalesourceov      = 'yes';
sourceZB                   = ft_sourceanalysis(cfg,zerobackdata);
sourceTB                   = ft_sourceanalysis(cfg, twobackdata);
sourceFX                   = ft_sourceanalysis(cfg, fixdata);

%% Visualize Source
bnd.pnt = sourcemodel.pos;
bnd.tri = sourcemodel.tri;
cfg = [];
cfg.projectmom = 'yes';

sourceDiff = sourceTB;
sourceDiff.avg.pow = sourceTB.avg.pow - sourceZB.avg.pow;

timepoint = [100 125 150 175 200 225];
tmpdatapows{1}=sourceZB.avg.pow;
tmpdatapows{2}=sourceTB.avg.pow;
tmpdatapows{3}=sourceFX.avg.pow;
tmpdatapows{4}=sourceDiff.avg.pow;
datalabels = {'Zero Back Powerspectrum', 'Two Back Powerspectrum', 'Fixation', 'Two Back - Zero Back Powerspectrum'};

for j =1:length(datalabels)
    figure;
    suptitle(datalabels{j});
    for i =1:length(timepoint)
        subplot(2,3, i);
        ft_plot_mesh(bnd, 'vertexcolor', tmpdatapows{j}(:,timepoint(i)));
        title([num2str(timepoint(i)),'ms']);
        colorbar;
    end
end


%% GMFP (Global mean field power)
% get channel mean
cfg                        = [];
cfg.trials                 = 1:length(data.trial);
cfg.channel                = 'MEG';
avgdata                    = ft_timelockanalysis(cfg, data);
datadiff = data;
for i=1:length(data.trial)
    datadiff.trial{i} = data.trial{i}-avgdata.avg;
end
