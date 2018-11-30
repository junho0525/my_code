function [source_ratio] = mnet_plot_band_contrast(config, fieldtripData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEEG Redefine Trials According to Band Power.                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.11.01 21:48 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input
% These five fields must be specified
if ~isfield(config, 'subjectID'),error('Invlid Subject ID!! There is no subjectID');end
if (~isfield(config, 'subjectPath')||~isdir(config.subjectPath));error('Invalid Subject Path!! There is no such folder');end
if (~isfield(config, 'savePath')||~isdir(config.savePath)),error('Invalid Save Path!! There is no such folder');end
if (~isfield(config, 'sessionID')),error('Invalid Session ID!! There is no session ID field!!');end
if (~isfield(config, 'band'))||~ischar(config.band)
    error('Invalid band!! There is no band field!!');
else
    config.band = lower(config.band);
    switch config.band
        case 'delta'
            config.bandfreq = [2,4];
        case 'theta'
            config.bandfreq = [5,7];
        case 'alpha'
            config.bandfreq = [8,12];
        case 'beta'
            config.bandfreq = [15,30];
        case 'gamma1'
            config.bandfreq = [30,60];
        case 'gamma2'
            config.bandfreq = [60,90];
        otherwise
            error('Invalid band!! You can specify band as ''delta'', ''theta'', ''alpha'', ''beta'', ''gamma1'', ''gamma2''');
    end
end
%% Set Path
rMEGPath = fullfile(config.subjectPath,[config.subjectID '_MEG_Restin_preproc'],config.subjectID,'MEG','Restin','rmegpreproc');
rMEGRawData=[config.subjectID '_MEG_' num2str(config.sessionID) '-Restin_rmegpreproc'];
anaPath = fullfile(config.subjectPath, [config.subjectID '_MEG_anatomy'],config.subjectID, 'MEG','anatomy');
headmodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_headmodel']);
sourcemodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_sourcemodel_2d.mat']);
% headmodelfile = 'E:\matlabworks\fieldtrip\fieldtrip-20180918\template\headmodel\standard_singleshell.mat';
% sourcemodelfile = 'E:\matlabworks\fieldtrip\fieldtrip-20180918\template\sourcemodel\cortex_8196.surf.gii';
%% Load MEEG Datas
use_my_library('ft',0);
use_my_library('spm',1);
if nargin==1 ||isempty(fieldtripData)
    error('There is no preprocessed fieldtrip file!!');
end
%% Prepare headmodel, sourcemodel, lead field
use_my_library('spm',0);
use_my_library('ft',1);
headmodel = load(headmodelfile);
if isfield(headmodel, 'headmodel')
    headmodel = headmodel.headmodel;
elseif isfield(headmodel, 'vol')
    headmodel = headmodel.vol;
end
[sourcePath sourceName ext]=fileparts(sourcemodelfile);
if strcmp(ext, '.mat')
    sourcemodel = load(sourcemodelfile);
    sourcemodel = sourcemodel.sourcemodel2d;
elseif strcmp(ext, '.gii')
    sourcemodel = ft_read_headshape(sourcemodelfile);
end
cfg = [];
cfg.grid = sourcemodel;
cfg.headmodel = headmodel;
cfg.channel = {'MEG'};
leadfield = ft_prepare_leadfield(cfg, fieldtripData.data);

%% Prepare Powerspectrum and find High / Low Alpha indices
%% Fourier analysis
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq = 1;
cfg.foilim = config.bandfreq; % set frequencies of interest as array
spectdata = ft_freqanalysis(cfg, fieldtripData.data);

%% Compute sensor level power spectra
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = config.bandfreq; % Alapha                     
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'yes';
datapow           = ft_freqanalysis(cfg, fieldtripData.data);
%% identify the indices of trials with high and low alpha power
freqind = nearest(datapow.freq, mean(config.bandfreq));
tmp     = datapow.powspctrm(:,:,freqind);
chanind = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));

clear tmp freqind
%% Check Power ratio 
%% compute the power spectrum for the median splitted data
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = (config.bandfreq(1)+config.bandfreq(2))/2;
cfg.trials       = indlow;
freq_low      = ft_freqanalysis(cfg, fieldtripData.data);

cfg.trials       = indhigh;
freq_high     = ft_freqanalysis(cfg, fieldtripData.data);

%% compute the beamformer filters based on the entire data
cfg                   = [];
cfg.frequency         = (config.bandfreq(1)+config.bandfreq(2))/2;
cfg.method            = 'pcc';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.fixedori      = 'yes';
source = ft_sourceanalysis(cfg, spectdata);

% use the precomputed filters
cfg                   = [];
cfg.frequency         = (config.bandfreq(1)+config.bandfreq(2))/2;
cfg.method            = 'pcc';
cfg.grid              = leadfield;
cfg.grid.filter       = source.avg.filter;
cfg.headmodel         = headmodel;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
source_low  = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_low));
source_high = ft_sourcedescriptives([], ft_sourceanalysis(cfg, freq_high));

cfg           = [];
cfg.operation = 'log10(x1)-log10(x2)';
cfg.parameter = 'pow';
source_ratio  = ft_math(cfg, source_high, source_low);
source_ratio.mask = (1+tanh(2.*(source_ratio.pow./max(source_ratio.pow(:))-0.5)))./2;
%% Source plot
cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = [-.3 .3];
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source_ratio);
% view([-90 30]);
light('style','infinite','position',[0 -200 200]);
title([config.subjectID '\_' num2str(config.sessionID)]);
