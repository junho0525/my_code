%% Load HCP resting state data, Time windows of interset
rmegpath = ('');
load([rmegpath '100307_MEG_3-Restin_rmegpreproc.mat']); 

%% Downsample
cfg = [];
cfg.resamplefs = 200;
cfg.detrend = 'yes';
datads = ft_resampledata(cfg, data);

%% use ICA in order to identify cardiac and blink components
cfg                 = [];
cfg.method          = 'runica'; 
cfg.runica.maxsteps = 50;
cfg.randomseed      = 0; % this can be uncommented to match the data that has been stored on disk
comp                = ft_componentanalysis(cfg, datads);

%% Reconstruct data from ica component
cfg = [];
cfg.component = [];% Specify badchannels;
dataica = ft_rejectcomponent(cfg, comp);

%% Calculate the Powersepctrum
cfg = [];
cfg.output = 'pow';
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq =1;
cfg.keeptrials = 'no';
datapow = ft_freqanalysis(cfg, dataica);


%% plot powerspectrum
figure;

% Delta (2-4Hz)
cfg        = [];
cfg.layout = '4D248_helmet.mat';
cfg.xlim   = [2 4];
subplot(3,2,1);
title('Delta (2-4Hz)');
ft_topoplotER(cfg, datapow);
% Theta (5-7Hz)
cfg        = [];
cfg.layout = '4D248_helmet.mat';
cfg.xlim   = [5 7];
subplot(3,2,2);
title('Theta (5-7Hz)');
ft_topoplotER(cfg, datapow);
% Alpha (8-12Hz)
cfg        = [];
cfg.layout = '4D248_helmet.mat';
cfg.xlim   = [8 12];
subplot(3,2,3);
title('Alpha (8-12Hz)');
ft_topoplotER(cfg, datapow);
% Beta (15-30Hz)
cfg        = [];
cfg.layout = '4D248_helmet.mat';
cfg.xlim   = [15 30];
subplot(3,2,4);
title('Beta (15-30Hz)');
ft_topoplotER(cfg, datapow);
% Gamma1 (30-60Hz)
cfg        = [];
cfg.layout = '4D248_helmet.mat';
cfg.xlim   = [30 60];
subplot(3,2,5);
title('Gamma1 (30-60Hz)');
ft_topoplotER(cfg, datapow);
% Gamma2 (60-90Hz)
cfg        = [];
cfg.layout = '4D248_helmet.mat';
cfg.xlim   = [60 90];
subplot(3,2,6);
title('Gamma2 (60-90Hz)');
ft_topoplotER(cfg, datapow);

%% Specify Forward MOdel. 
headmodel = load([rmegpath '100307_MEG_anatomy_headmodel']);
headmodel = headmodel.headmodel;
sourcemodel = load([rmegpath '100307_MEG_anatomy_sourcemodel_2d']);
sourcemodel = sourcemodel.sourcemodel2d;

% % visualize the coregisteration of sensor, headmodel, and sourcemodel.
% figure;
% 
% % make the headmodel surface transparent
% ft_plot_vol(hdm, 'edgecolor', 'none'); alpha 0.4           
% ft_plot_mesh(ft_convert_units(sourcemodel, 'cm'),'vertexcolor',sourcemodel.sulc);
% ft_plot_sens(dataclean.grad);
% view([0 -90 0])

% Calculate Lead Field
cfg = [];
cfg.grid = sourcemodel;
cfg.headmodel = headmodel;
cfg.channel = {'MEG'};
leadfield = ft_prepare_leadfield(cfg, dataica);

%% Frequency Analysis for Source Reconstruction
cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq = 1;
cfg.foi = 10; % set frequencies of interest as array
freq = ft_freqanalysis(cfg, dataica);

%% Source Reconstruction
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.fixedori      = 'yes';
source = ft_sourceanalysis(cfg, freq);
source = ft_sourcedescriptives([], source); % to get the neural-activity-index

%% plot the neural activity index (power/noise)
cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'nai';
cfg.maskparameter = cfg.funparameter;
cfg.funcolorlim   = 'auto'; % [0.0 8]
cfg.opacitylim    = [3 8]; 
cfg.opacitymap    = 'rampup';  
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source);
view([-90 30]);
light;

%% compute sensor level single trial power spectra
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = [8 12];                          
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'yes';
datapow           = ft_freqanalysis(cfg, dataica);

%% identify the indices of trials with high and low alpha power
freqind = nearest(datapow.freq, 10);
tmp     = datapow.powspctrm(:,:,freqind);    
chanind = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));
%% compute the power spectrum for the median splitted data
cfg              = [];
cfg.trials       = indlow; 
datapow_low      = ft_freqdescriptives(cfg, datapow);

cfg.trials       = indhigh; 
datapow_high     = ft_freqdescriptives(cfg, datapow);

%% compute the difference between high and low
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'divide';
powratio      = ft_math(cfg, datapow_high, datapow_low);

%% plot the topography of the difference along with the spectra
cfg        = [];
cfg.layout = '4D248_helmet.mat';
cfg.xlim   = [9.9 10.1];
figure; ft_topoplotER(cfg, powratio);

cfg         = [];
cfg.channel = {'A196'};
figure; ft_singleplotER(cfg, datapow_high, datapow_low);

%% compute fourier spectra for frequency of interest according to the trial split
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'fourier';
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 1;
cfg.foi        = 10;

cfg.trials = indlow; 
freq_low   = ft_freqanalysis(cfg, dataica);

cfg.trials = indhigh; 
freq_high  = ft_freqanalysis(cfg, dataica);

%% compute the beamformer filters based on the entire data
cfg                   = [];
cfg.frequency         = freq.freq;
cfg.method            = 'pcc';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.keeptrials        = 'yes';
cfg.pcc.lambda        = '10%';
cfg.pcc.projectnoise  = 'yes';
cfg.pcc.keepfilter    = 'yes';
cfg.pcc.fixedori      = 'yes';
source = ft_sourceanalysis(cfg, freq);

% use the precomputed filters 
cfg                   = [];
cfg.frequency         = freq.freq;
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

% create a fancy mask
source_ratio.mask = (1+tanh(2.*(source_ratio.pow./max(source_ratio.pow(:))-0.5)))./2; 

cfg = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = 'mask';
cfg.funcolorlim   = [-.3 .3];
cfg.funcolormap   = 'jet';
cfg.colorbar      = 'no';
ft_sourceplot(cfg, source_ratio);
view([-90 30]);
light('style','infinite','position',[0 -200 200]);
