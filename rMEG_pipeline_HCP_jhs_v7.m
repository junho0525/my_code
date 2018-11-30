function [D,DCM] = rMEG_pipeline_HCP_jhs_v7(config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCP Rest MEG CSD DCM Pipeline Using SPM                                 %
%                                                                         %
% This pipeline code is for rest HCP dataset                              %
% This pipeline split megdata into high/low alpha                         %
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
%    DCM parameters                                                       %
%       config.TimeWin    - must need
%       config.analysis   - default('CSD')
%       config.ROI       
%       config.ROIName    
%       config.effect     - optional
%       config.effectName - optional
%       config.trials     - set to 1 for continuous data
%       config.ChanMode     - number of channel modes
%       config.A
%       config.B
%       config.C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.05.03 00:31 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input
% These five fields must be specified
if ~isfield(config, 'subjectID'),error('Invlid Subject ID!! There is no subjectID');end
if ~isfield(config, 'taskType'),error('Invalid Task Type!! There is no taskType');end
if (~isfield(config, 'subjectPath')||~isdir(config.subjectPath));error('Invalid Subject Path!! There is no such folder');end
if (~isfield(config, 'savePath')||~isdir(config.savePath)),error('Invalid Save Path!! There is no such folder');end
if (~isfield(config, 'sessionID')),error('Invalid Session ID!! There is no session ID field!!');end

% DCM parameters
if ~isfield(config, 'TimeWin'), error('you need specify time window to analyse'); end % default : 0
if ~isfield(config, 'FreqWin'), config.FreqWin = [1, 100]; end % default : [1,100]
if ~isfield(config, 'analysis'), config.analysis = 'CSD'; end % default : 'CSD'
if ~isfield(config, 'model'), config.model = 'ERP'; end % default : 'EPR'
if ~isfield(config, 'spatial'), config.spatial = 'ECD'; end % default : 'ECD'
if ~isfield(config, 'ROI')
    config.ROI     = [[-46;-66;30] [49;-63;33] [0; -58; 0] [-1;54;27]];
    config.ROIName = {'lLP', 'rLP', 'Prec', 'mPFC'};
end
if ~isfield(config, 'trials'), config.trials =1; end
if ~isfield(config, 'chanMode'), config.chanMode = 8; end % default : 8
if ~isfield(config, 'A')
    config.A = 'default'; % A{1} : Forward, A{2} : Backward, A{3} : Modulatory
end
if ~isfield(config, 'B')
    config.B = 'default'; % B : Effect dependent modulation
end
if ~isfield(config, 'C')  
    switch config.analysis
        case 'CSD' % 'CSD' do not have any input
        otherwise
            config.C = 'default'; % C : input vector
    end
end
if ~isfield(config, 'effect')
    config.effect = 'default';
    config.effectName = {};
else
    if ~isfield(config, 'effectName')
        config.effectName = cell(1, size(config.effect,2));
        for i = 1:size(config.effect,2)
            config.effectName{i} = ['effect ' num2str(i)]; 
        end
    elseif length(config.effectName)~=size(config.effect,2)
        config.effectName = [];
        for i = 1:size(config.effect,2)
            config.effectName{i} = ['effect ' num2str(i)]; 
        end
    end
end

% These parameters are optional
if ~isfield(config, 'showMEGSignal'), config.showMEGSignal = 0; end
if ~isfield(config, 'showForwardResult'), config.showForwardResult =0; end
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
try
    D = spm_eeg_load(fullfile(config.savePath,['afdfspm_', rMEGRawData]));
    fieldtripData = load(fullfile(rMEGPath,rMEGRawData));
catch
    fieldtripData = load(fullfile(rMEGPath, rMEGRawData));
    try
        D = spm_eeg_load(config.savePath, ['spm_',rMEGRawData]);
    catch
        D = spm_eeg_ft2spm(fieldtripData.data, fullfile(config.savePath,['spm_' ,rMEGRawData]));
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
%% Prepare Powerspectrum and find High / Low Alpha indices
%% Compute sensor level power spectra
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.taper        = 'dpss';
cfg.foilim       = [8 12]; % Alapha                     
cfg.tapsmofrq    = 1;             
cfg.keeptrials   = 'yes';
datapow           = ft_freqanalysis(cfg, fieldtripData.data);
%% identify the indices of trials with high and low alpha power
freqind = nearest(datapow.freq, 10);
tmp     = datapow.powspctrm(:,:,freqind);    
chanind = find(mean(tmp,1)==max(mean(tmp,1)));  % find the sensor where power is max
indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));
clear tmp freqind chanind
%% Redefine trials
condlength = length(D.conditions);
conditionlabels = cell(1,condlength);
for i = 1:147
    if ~isempty(find(~(indlow-i), 1))
        conditionlabels{i} = 'Alpha Low';
    elseif ~isempty(find(~(indhigh-i), 1))
        conditionlabels{i} = 'Alpha High';
    end
end
D = conditions(D,':',conditionlabels);
save(fullfile(config.savePath,D.fname),'D');
%% Check MEG signals
if config.showMEGSignal
    figure;
    suptitle('All MEG Signals use trial 1');
    for i = 1:length(D.condlist)
        subplot(length(D.condlist),1,i);
        conditionIdx = find(strcmp(D.conditions,D.condlist{i}));
        for j = conditionIdx
            plot(D.time,D(1,:,j));
            hold on;
            for k = 2:size(D,1)
                plot(D.time,D(k,:,j));
            end
            hold off;
        end
        title(['MEG Signal - condition ' D.condlist{i}]);
    end
end
%%
%% Check forward
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
    load(sourcemodelfile); % This file is originally in anatomy directory
    sourcemodel2d=ft_convert_units(sourcemodel2d, 'mm');
    sourcemodelsubj = sourcemodel2d;

    load(headmodelfile);
    headmodel = ft_convert_units(headmodel, 'cm');

    grad=fieldtripData.data.grad;

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
    job.vol = headmodel;
    job.grid = sourcemodelsubj;
    job.grad = grad;
    job.channel = channels;
    job.normalize = 'yes';
    job.reducerank = 2;

    val = []; tlck=[]; sourcemodelsubj = []; resultprefix = []; options = []; mixing = []; megpath = []; rmeg = []; mainpath = [];
    i = []; headmodel = []; grad = []; ft_default = []; experimentid = []; comp_class = []; channels = []; anaPath = [];
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
%% DCM %%
%% Specify parameters
%% Prepare data ( Set DCM.xY )
DCM = [];
DCM.xY.Dfile = D.fullfile;
DCM.xY.modality = 'MEG';
DCM.name = ['DCM_' D.fname];
DCM.val =D.val;
DCM.options.analysis = config.analysis;
DCM.options.model = config.model;
DCM.options.spatial = config.spatial;
DCM.options.trials = config.trials;
DCM.options.Tdcm = config.TimeWin;
DCM.options.Fdcm = config.FreqWin;
DCM.options.Nmodes = config.chanMode;
DCM.options.h = 1;
DCM.options.D = 1;
DCM.options.lock = 1;
DCM.options.multiC = 0;
DCM.options.symmetry = 0;
DCM.options.Rft = 5;
DCM = spm_dcm_erp_data(DCM, 0);
spm_dcm_erp_results(DCM, 'Data');

%% Prepare dipolefit. ( Set DCM.M . Only need for 'ECD' and 'IMG'. No need for 'LFP)
DCM.Lpos   = config.ROI;
DCM.Sname  = config.ROIName;
% DCM.Lpos = [ [-46;19;22] [41;31;30] [-31;21;4] [34;23;1]];
% DCM.Sname = {'lDLPFC','rDLPRF','lVLPRF','rVLPFC'};
Nareas = size(DCM.Lpos,2);
% Spatial Model
DCM = spm_dcm_erp_dipfit_jhs(DCM);

%% Prepare Connectivity Matrices ( Set A, B, C - No need C for 'CSD' )
if ischar(config.effect)
    DCM.xU.name = {}; % If there is some between trial effect, then you should specify xU.
    DCM.xU.X    = zeros(length(D.condlist),0); % Design Matrix
else
    if size(config.effect,1) == length(D.condlist)
        DCM.xU.name = config.effectName;
        DCM.xU.X    = config.effect;
    else
        warning('Size of Design Matrix does not match with the number of conditions');
        warning('By Default, we do not use Design matrix.');
        DCM.xU.name = {}; % If there is some between trial effect, then you should specify xU.
        DCM.xU.X    = zeros(length(D.condlist),0); % Design Matrix
    end
end




% Specify connectivity model ( Set as full model )
if ischar(config.A)
    DCM.A{1} = ones(Nareas, Nareas); % A{1} : Forward
    DCM.A{2} = ones(Nareas,Nareas); % A{2} : Backward
    DCM.A{3} = zeros(Nareas,Nareas); % A{3} : Modulatory
elseif iscell(config.A)
    flag = [];
    for i = 1:3
        flag = [flag find(size(config.A{i})~=[Nareas, Nareas])];
    end
    if ~isempty(flag)
        error('Invalid A marix size!!');
    end
    DCM.A = config.A;
end
if ischar(config.B)
    DCM.B{1} = zeros(Nareas,Nareas); % B{1} : Effect dependent modulatory
elseif iscell(config.B)
    flag = [];
    for i = 1:length(config.B)
        flag = [flag find(size(config.A{i})~=[Nareas, Nareas])];
    end
    if ~isempty(flag)
        error('Invalid B marix size!!');
    end
    DCM.B = config.B;
end
% C matrix is not used in CSD.
%% Invert
DCM = spm_dcm_csd_jhs(DCM);
