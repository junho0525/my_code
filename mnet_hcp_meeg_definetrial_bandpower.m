function [D,fieldtripData] = mnet_hcp_meeg_definetrial_bandpower(config, D, fieldtripData)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEEG Redefine Trials According to Band Power.                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.05.14 12:12 - By Junho Son                                     %
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
sourcemodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_sourcemodel_2d']);
%% Load MEEG Datas
use_my_library('ft',0);
use_my_library('spm',1);
if nargin==1 || isempty(D)||isempty(fieldtripData)
    try
        D = spm_eeg_load(fullfile(config.savePath,['afdfspm_', rMEGRawData]));
        fieldtripData = load(fullfile(rMEGPath,rMEGRawData));
    catch
        error(['There is no preprocessed file!! filename : ' fullfile(config.savePath,['afdfspm_', rMEGRawData])]);
    end
end
% if length(D.condlist)>1
%     condpart = strsplit(D.condlist{1}, ' ');
%     if strcmp(lower(condpart{1}),config.band)
%         return; % If already calculated, then return immediately.
%     end
% end
%% Prepare headmodel, sourcemodel, lead field
use_my_library('spm',0);
use_my_library('ft',1);
headmodel = load(headmodelfile);
headmodel = headmodel.headmodel;
sourcemodel = load(sourcemodelfile);
sourcemodel = sourcemodel.sourcemodel2d;

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

%% Source Reconsturction
cfg                   = [];
cfg.frequency         = config.bandfreq;
cfg.method            = 'dics';
cfg.grid              = leadfield;
cfg.headmodel         = headmodel;
cfg.rawtrial          = 'yes';
cfg.dics.lambda       = 0;
cfg.dics.projectnoise  = 'yes';
source = ft_sourceanalysis(cfg, spectdata);
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

for i =1:length(datapow.freq)
    tmp     = tmp + datapow.powspctrm(:,:,i);
end
chanind = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
% indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
% indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));

% split trials using ROI power
use_my_library('ft',0);
use_my_library('spm',1);
tempDistance = D.inv{1}.mesh.tess_mni.vert - repmat([0,-58,0],8004,1);
tempDistance = sum(tempDistance.^2,2);
roiInd = find(min(tempDistance)==tempDistance);
tmp = [];
for i=1:length(source.trial)
    tmp(end+1) = source.trial(i).pow(roiInd);
end
indlow  = find(tmp(:)<=median(tmp(:)));
indhigh = find(tmp(:)>=median(tmp(:)));

cfg              = [];
cfg.trials       = indlow; 
datapow_low      = ft_freqdescriptives(cfg, datapow);

cfg.trials       = indhigh; 
datapow_high     = ft_freqdescriptives(cfg, datapow);

clear tmp freqind
%% Redefine trials
condlength = length(D.conditions);
conditionlabels = cell(1,condlength);
for i = 1:147
    if ~isempty(find(~(indlow-i), 1))
        conditionlabels{i} = [upper(config.band(1)) config.band(2:end) ' Low'];
    elseif ~isempty(find(~(indhigh-i), 1))
        conditionlabels{i} = [upper(config.band(1)) config.band(2:end) ' High'];
    end
end
D = conditions(D,':',conditionlabels);
save(fullfile(config.savePath,D.fname),'D');
%% Check Power ratio 
if config.showPowerRatio
    use_my_library('spm',0);
    use_my_library('ft',1);
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
    cfg.frequency         = (config.bandfreq(1)+config.bandfreq(2))/2;;
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
%     %% compute the difference between high and low
%     cfg = [];
%     cfg.parameter = 'powspctrm';
%     cfg.operation = 'divide';
%     powratio      = ft_math(cfg, datapow_high, datapow_low);
    %% Source plot
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
    %% plot the topography of the difference along with the spectra
%     cfg        = [];
%     cfg.layout = '4D248_helmet.mat';
%     cfg.xlim   = config.bandfreq;
%     figure; ft_topoplotER(cfg, powratio);

    cfg         = [];
    cfg.channel = {datapow.label{chanind}};
    figure; ft_singleplotER(cfg, datapow_high, datapow_low);
    use_my_library('ft',0);
    use_my_library('spm',1);
    clear chanind
end