function [D,DCM] = power_ratio_map_one_source(config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCP Rest MEG CSD DCM Pipeline Using SPM                                 %
%                                                                         %
% This pipeline code is for rest HCP dataset                              %
% This pipeline split megdata into source level high/low frequency band   %
%                                                                         %
% [D, DCM] = rMEG_pipeline_HCP_jhs_v6(config)                             %
%    Output      D     - spm meeg object                                  %
%                DCM   - spm DCM struct                                   %
%                                                                         %
%    Input                                                                %
%    General parameters                                                   %
%       config.subjectID     - HCP dataset subject id                     %
%       config.sessionID     - HCP dataset session id                     %
%       config.tasktype      - 'Wrkmem', 'Motort', 'StoryM'               %
%       config.subjectPath   - HCP Single Subject Data Path               %
%       config.savePath      - Save path for output and intermediate files%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.06.04 11:03 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input
% These five fields must be specified
if ~isfield(config, 'subjectID'),error('Invlid Subject ID!! There is no subjectID');end
if ~isfield(config, 'taskType'),error('Invalid Task Type!! There is no taskType');end
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

% These parameters are optional
if ~isfield(config, 'showForwardResult'), config.showForwardResult =0; end
if ~isfield(config, 'showPowerRatio'),config.showPowerRatio =0; end

%% Set Path
switch config.taskType
    case 'Restin'
        rMEGPath = fullfile(config.subjectPath,[config.subjectID '_MEG_Restin_preproc'],config.subjectID,'MEG','Restin','rmegpreproc');
        rMEGRawData=[config.subjectID '_MEG_' num2str(config.sessionID) '-Restin_rmegpreproc'];
        icaInfo = fullfile(config.subjectPath,[config.subjectID '_MEG_Restin_preproc'],config.subjectID,'MEG','Restin','icaclass',[config.subjectID,'_MEG_' num2str(config.sessionID) '-Restin_icaclass_vs.mat']);
    otherwise
        Error();
end
anaPath = fullfile(config.subjectPath, [config.subjectID '_MEG_anatomy'],config.subjectID, 'MEG','anatomy');
headmodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_headmodel']);
sourcemodelfile = fullfile(anaPath,[config.subjectID '_MEG_anatomy_sourcemodel_2d']);
%% Load fieldtrip data
use_my_library('ft',0);
use_my_library('spm',1);
try
    D = spm_eeg_load(fullfile(config.savePath,['afdfspm_', rMEGRawData]));
    fieldtripData = load(fullfile(rMEGPath,rMEGRawData));
    if isfield(fieldtripData, 'data'), fieldtripData = fieldtripData.data; end
catch
    fieldtripData = load(fullfile(rMEGPath, rMEGRawData));
    if isfield(fieldtripData, 'data'), fieldtripData = fieldtripData.data; end
    try
        D = spm_eeg_load(config.savePath, ['spm_',rMEGRawData]);
    catch
        D = spm_eeg_ft2spm(fieldtripData, fullfile(config.savePath,['spm_' ,rMEGRawData]));
    end
    %% 1-20Hz bandpass filter & downsampling
    S.D = D;
    S.band = 'high';
    S.freq = 0.1;
    D = spm_eeg_filter(S);
    S = [];
    S.D = D;
    S.method = 'fft';
    S.fsample_new = 400;
    D = spm_eeg_downsample(S);
    S = [];
    S.D = D;
    S.band = 'low';
    S.freq = 200;
    D = spm_eeg_filter(S);
    %% Artefact detaction
    S = [];
    S.D = D;
    S.mode = 'reject';
    S.methods.channels = {'MEG'};
    S.methods.fun = 'threshchan';
    S.methods.settings.threshold = 3000;%(3000 fT)
    S.methods.settings.excwin = 1000;
    D = spm_eeg_artefact(S);
    delete(fullfile(config.savePath,['fdfspm_', rMEGRawData, '.mat']));
    delete(fullfile(config.savePath,['fdfspm_', rMEGRawData, '.dat']));
    delete(fullfile(config.savePath,['dfspm_', rMEGRawData, '.mat']));
    delete(fullfile(config.savePath,['dfspm_', rMEGRawData, '.dat']));
    delete(fullfile(config.savePath,['fspm_', rMEGRawData, '.mat']));
    delete(fullfile(config.savePath,['fspm_', rMEGRawData, '.dat']));
    delete(fullfile(config.savePath,['spm_', rMEGRawData, '.mat']));
    delete(fullfile(config.savePath,['spm_', rMEGRawData, '.dat']));
end
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
leadfield = ft_prepare_leadfield(cfg, fieldtripData);

%% Check and specify Forward Model
use_my_library('ft',0);
use_my_library('spm',1);
if isfield(D,'inv') &&isfield(D.inv{1},'forward')&&isfield(D.inv{1},'datareg')&&isfield(D.inv{1},'mesh')
else
    %% Specify Forward Model
    load(icaInfo); % This file is originally in Restin/icaclass directory of Restin_preproc
    % read the source and volume conduction model from current dir with
    % outputs of previous pipelines
    fid = fopen(fullfile(anaPath,[config.subjectID '_MEG_anatomy_transform.txt']), 'rt');
    strFid=fread(fid,[1 inf], '*char');
    eval(strFid);
    fclose(fid);
    fid = [];
    clear fid;
    strFid = [];
    clear strFid;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % specify sourcemodel, headmodel, sensor and transform if it necessary
    sourcemodelsubj=ft_convert_units(sourcemodel, 'mm');
    headmodelsubj = ft_convert_units(headmodel, 'cm');

    grad=fieldtripData.grad;

    % These three code lines are for EEG
    % gradBalanced = grad;
    % gradBalanced = ft_apply_montage(gradBalanced, gradBalanced.balance.Supine, 'keepunused', 'yes', 'inverse', 'yes');
    % grad=gradBalanced;

    grad = ft_convert_units(grad, 'cm');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare the component data in order for ft_sourceanalysis to be able to
    % swallow it
    channels = comp_class.topolabel;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the forward solution

    job = [];
    job.vol = headmodelsubj;
    job.grid = sourcemodelsubj;
    job.grad = grad;
    job.channel = channels;
    job.normalize = 'yes';
    job.reducerank = 2;

    val = []; tlck=[]; sourcemodelsubj = []; resultprefix = []; options = []; mixing = []; megpath = []; rmeg = []; mainpath = [];
    i = []; headmodelsubj = []; grad = []; ft_default = []; experimentid = []; comp_class = []; channels = []; anaPath = [];
    clear val tlck sourcemodelsubj resultprefix options mixing megpath;
    clear rmeg mainpath i headmodel grad ft_default experimentid comp_class channels anaPath;

    % specify forward model (Bring forward model information 
    D.inv = [];
    D.inv{1}.forward = [];
    D.inv{1}.method='Imaging';
    D.inv{1}.forward.toMNI = transform.bti2spm;
    D.inv{1}.forward.fromMNI = transform.spm2bti;
    D.inv{1}.forward.voltype='Single Shell';
    D.inv{1}.forward.vol = ft_convert_units(job.vol,'m');
    D.inv{1}.forward.modality = 'MEG';
    D.inv{1}.forward.siunits=1;
    D.inv{1}.forward.mesh=[];
    D.inv{1}.forward.mesh_correction=[];
    D.inv{1}.forward.mesh.face = job.grid.tri;
    grid = ft_convert_units(job.grid,'m');
    D.inv{1}.forward.mesh.vert = grid.pos;
    D.inv{1}.forward.sensors = [];
    D.inv{1}.forward.sensors = ft_convert_units(job.grad,'m');
    D.inv{1}.mesh.tess_mni.face = D.inv{1}.forward.mesh.face;
    grid = ft_convert_units(job.grid,'mm');
    D.inv{1}.mesh.tess_mni.vert = spm_eeg_inv_transform_points(D.inv{1}.forward.toMNI,grid.pos); % tess_mni should have mni coordinate
    D.inv{1}.forward.toMNI = [D.inv{1}.forward.toMNI(:, 1:3)*1000 D.inv{1}.forward.toMNI(:,4)];
    D.inv{1}.forward.fromMNI = eye(4)/D.inv{1}.forward.toMNI
    D.inv{1}.datareg.modality = 'MEG';
    D.inv{1}.datareg.sensors = ft_convert_units(job.grad,'m');
    D.inv{1}.datareg.toMNI = D.inv{1}.forward.toMNI;
    D.inv{1}.datareg.fromMNI = D.inv{1}.forward.fromMNI;
    Datareg = []; grid = [];
    save(fullfile(config.savePath,D.fname),'D');
    clear Datareg grid;
    % check meshes and co-registration
    if config.showForwardResult; spm_eeg_inv_checkforward(D); end
    save(fullfile(config.savePath,D.fname),'D');
end
%% Prepare Powerspectrum and find High / Low Alpha indices
%% Fourier analysis
use_my_library('spm',0);
use_my_library('ft',1);

cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = config.bandfreq;
cfg.tapsmofrq    = 1;
cfg.keeptrials   = 'yes';
datapow           = ft_freqanalysis(cfg, fieldtripData);

use_my_library('ft',0);
use_my_library('spm',1);
%% identify the indices of trials with high and low alpha power
freqind = nearest(datapow.freq, 10);
tmp     = datapow.powspctrm(:,:,freqind);    
chanind = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));


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
    freq_low      = ft_freqanalysis(cfg, fieldtripData);

    cfg.trials       = indhigh; 
    freq_high     = ft_freqanalysis(cfg, fieldtripData);

    %% Fourier analysis
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.output = 'fourier';
    cfg.keeptrials = 'yes';
    cfg.tapsmofrq = 1;
    cfg.foilim = config.bandfreq; % set frequencies of interest as array
    spectdata = ft_freqanalysis(cfg, fieldtripData);
    %% compute the beamformer filters based on the entire data
    cfg                   = [];
    cfg.frequency         = config.bandfreq;
    cfg.method            = 'dics';
    cfg.grid              = leadfield;
    cfg.headmodel         = headmodel;
    cfg.dics.lambda       = '10%';
    cfg.dics.projectnoise = 'yes';
    cfg.dics.keepfilter   = 'yes';
    source = ft_sourceanalysis(cfg, spectdata);

    % use the precomputed filters 
    cfg                   = [];
    cfg.frequency         = config.bandfreq;
    cfg.method            = 'dics';
    cfg.grid              = leadfield;
    cfg.grid.filter       = source.avg.filter;
    cfg.headmodel         = headmodel;
    cfg.dics.lambda       = '10%';
    cfg.dics.projectnoise = 'yes';
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
    view([-90 30]);
    light('style','infinite','position',[0 -200 200]);
    %% Compute sensor level power spectra
    cfg               = [];
    cfg.output        = 'pow';
    cfg.method        = 'mtmfft';
    cfg.taper         = 'dpss';
    cfg.foilim        = config.bandfreq; % Alapha                     
    cfg.tapsmofrq     = 1;             
    cfg.keeptrials    = 'yes';
    datapow           = ft_freqanalysis(cfg, fieldtripData);
    cfg               = [];
    cfg.trials        = indlow; 
    datapow_low       = ft_freqdescriptives(cfg, datapow);

    cfg.trials        = indhigh; 
    datapow_high      = ft_freqdescriptives(cfg, datapow);
%     cfg        = [];
%     cfg.layout = '4D248_helmet.mat';
%     cfg.xlim   = config.bandfreq;
%     figure; ft_topoplotER(cfg, powratio);

    figure;
    plot(datapow.freq, mean(datapow_high.powspctrm,1));
    hold on
    plot(datapow.freq, mean(datapow_low.powspctrm,1));
    hold off
    use_my_library('ft',0);
    use_my_library('spm',1);
    clear chanind
end