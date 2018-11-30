%% Load HCP resting state data, Time windows of interset
load('100307_MEG_3-Restin_rmegpreproc.mat');
D = spm_eeg_load('fdf_spm_100307_MEG_3-Restin_rmegpreproc');
cfg.toilim = [-0.5 1.5];
data = ft_redefinetrial(cfg, data); % Time windows of interest
%% Calculating the Cross Spectral Density Matrix
cfg = [];
cfg.method = 'mtmfft';
cfg.output    = 'powandcsd';
cfg.tapsmofrq = 4;
cfg.foilim    = [2 60];
freqData = ft_freqanalysis(cfg,data);

%% Compute LeadField
dippos = D.inv{1}.forward.mesh.vert;
hcp_read_matlab('100307_MEG_anatomy_headmodel.mat');
grad = D.inv{1}.forward.sensors;
headmodel = ft_convert_units(headmodel , 'm');
hcp_read_matlab('100307_MEG_anatomy_sourcemodel_2d');
cfg.grad = grad;
cfg.grid = sourcemodel2d;
cfg.headmodel = headmodel;
cfg.channel = D.chanlabels;
grid = ft_prepare_leadfield(cfg);

%% Invert Solution

cfg.grid = grid;
cfg.method = 'dics';
cfg.frequency = 18;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda = 0;
source = ft_sourceanalysis(cfg,data);