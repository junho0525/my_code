function [D,DCM_low, DCM_high] = mnet_HCPrMEG_AlphaContrastDCM(cfg, D, fieldtripData, sourcemodelfile, headmodelfile,transformfile)
if ~isfield(cfg.definetrial, 'chanlabel')
    cfg.definetrial.chanlabel = 'max';
end

%% Preprocess
% Bandpass filter & downsampling
rMEGRawData = strsplit(D.fname,'.');
rMEGRawData = rMEGRawData{1};

S = [];
S.D = D;
S.method = 'resample';
S.fsample_new = cfg.preprocess.downsample;
if cfg.preprocess.downsample~=D.fsample
    D = spm_eeg_downsample(S); % resample of time has some edge effect problem
end

S = [];
S.D = D;
S.band = 'high';
S.freq = cfg.preprocess.freqband(1);
D = spm_eeg_filter(S);

S = [];
S.D = D;
S.band = 'low';
S.freq = cfg.preprocess.freqband(2);
D = spm_eeg_filter(S);
% Artefact detaction
S = [];
S.D = D;
S.mode = 'reject';
S.methods.channels = {'MEG'};
S.methods.fun = 'threshchan';
S.methods.settings.threshold = 3000;%(3000 fT)
S.methods.settings.excwin = 1000;
D = spm_eeg_artefact(S);

delete(fullfile(cfg.savePath,['ffd',rMEGRawData, '.mat']));
delete(fullfile(cfg.savePath,['ffd', rMEGRawData, '.dat']));
delete(fullfile(cfg.savePath,['fd', rMEGRawData, '.mat']));
delete(fullfile(cfg.savePath,['fd', rMEGRawData, '.dat']));
delete(fullfile(cfg.savePath,['d', rMEGRawData, '.mat']));
delete(fullfile(cfg.savePath,['d', rMEGRawData, '.dat']));
delete(fullfile(cfg.savePath,[rMEGRawData, '.mat']));
delete(fullfile(cfg.savePath,[rMEGRawData, '.dat']));


%% Set Forward Model
if isfield(D,'inv') &&isfield(D.inv{1},'forward')&&isfield(D.inv{1},'datareg')&&isfield(D.inv{1},'mesh')
else
    % Specify Forward Model
    % read the source and volume conduction model from current dir with
    % outputs of previous pipelines
    fid = fopen(transformfile, 'rt');
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
    channels = D.chanlabels;

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
    save(fullfile(cfg.savePath,D.fname),'D');
    clear Datareg grid;
    % check meshes and co-registration
    % if cfg.showForwardResult; spm_eeg_inv_checkforward(D); end
    save(fullfile(cfg.savePath,D.fname),'D');
end
%% Redefine Trials As High / Low Alpha
% Compute sensor level power spectra
use_my_library('spm',0);
use_my_library('ft',1);

tmpcfg              = [];
tmpcfg.output       = 'pow';
tmpcfg.method       = 'mtmfft';
tmpcfg.taper        = 'dpss';
tmpcfg.foilim       = cfg.definetrial.bandfreq; % Alapha                     
tmpcfg.tapsmofrq    = 1;             
tmpcfg.keeptrials   = 'yes';
datapow           = ft_freqanalysis(tmpcfg, fieldtripData.data);
% identify the indices of trials with high and low alpha power
freqind = nearest(datapow.freq, mean(cfg.definetrial.bandfreq));
tmp     = datapow.powspctrm(:,:,freqind);
if strcmp(cfg.definetrial.chanlabel, 'max') || isempty(find(ismember(D.chanlabels,cfg.definetrial.chanlabel))) % find the sensor where power is max (Default)
    chanind = find(mean(tmp,1)==max(mean(tmp,1)));  
    indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
    indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));
else
    chanind = find(ismember(D.chanlabels,cfg.definetrial.chanlabel));
    indlow  = find(tmp(:,chanind)<=median(tmp(:,chanind)));
    indhigh = find(tmp(:,chanind)>=median(tmp(:,chanind)));
end
tmpcfg              = []; 
tmpcfg.trials       = indlow; 
datapow_low      = ft_freqdescriptives(tmpcfg, datapow);

tmpcfg.trials       = indhigh; 
datapow_high     = ft_freqdescriptives(tmpcfg, datapow);

clear tmp freqind
% Redefine trials
use_my_library('ft',0);
use_my_library('spm',1);
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
save(fullfile(cfg.savePath,D.fname),'D');

%% Run spDCM
%  show ROIs
dip.n_seeds = 1;
dip.n_dip = size(cfg.DCM.ROI,2);
dip.Mtb=1;
dip.j{1} = zeros(3*dip.n_dip,1);
dip.loc{1} = cfg.DCM.ROI;
spm_eeg_inv_ecd_DrawDip('init',dip);

%% DCM Low Alpha
DCM_low = [];
DCM_low.xY.Dfile = D.fullfile;
DCM_low.xY.modality = 'MEG';
DCM_low.name = ['DCM_' D.fname '_LowAlpha'];
DCM_low.val =D.val;
DCM_low.options.analysis = cfg.DCM.analysis;
DCM_low.options.model = cfg.DCM.model;
DCM_low.options.spatial = cfg.DCM.spatial;
DCM_low.options.trials = cfg.DCM.trials{1};
DCM_low.options.Tdcm = cfg.DCM.TimeWin;
DCM_low.options.Fdcm = cfg.DCM.FreqWin;
DCM_low.options.Nmodes = cfg.DCM.chanMode;
DCM_low.options.h = 1;
DCM_low.options.D = 1;
DCM_low.options.lock = 1;
DCM_low.options.multiC = 0;
DCM_low.options.symmetry = 0;
DCM_low.options.Rft = 5;
DCM_low = spm_dcm_erp_data(DCM_low, 0); % 0 for not to make erp
%spm_dcm_erp_results(DCM, 'Data');

% Prepare dipolefit. ( Set DCM.M . Only need for 'ECD' and 'IMG'. No need for 'LFP)
DCM_low.Lpos   = cfg.DCM.ROI;
DCM_low.Sname  = cfg.DCM.ROIName;
Nareas = size(DCM_low.Lpos,2);
% Spatial Model
DCM_low = spm_dcm_erp_dipfit_jhs(DCM_low);

% Prepare Connectivity Matrices ( Set A, B, C - No need C for 'CSD' )
if ischar(cfg.DCM.effect)
    DCM_low.xU.name = {}; % If there is some between trial effect, then you should specify xU.
    DCM_low.xU.X    = zeros(length(D.condlist),0); % Design Matrix
else
    if size(cfg.DCM.effect,1) == length(D.condlist)
        DCM_low.xU.name = cfg.DCM.effectName;
        DCM_low.xU.X    = cfg.DCM.effect;
    else
        warning('Size of Design Matrix does not match with the number of conditions');
        warning('By Default, we do not use Design matrix.');
        DCM_low.xU.name = {}; % If there is some between trial effect, then you should specify xU.
        DCM_low.xU.X    = zeros(length(D.condlist),0); % Design Matrix
    end
end
% Specify connectivity model ( Set as full model )
if ischar(cfg.DCM.A)
    DCM_low.A{1} = ones(Nareas, Nareas); % A{1} : Forward
    DCM_low.A{2} = ones(Nareas,Nareas); % A{2} : Backward
    DCM_low.A{3} = zeros(Nareas,Nareas); % A{3} : Modulatory
elseif iscell(cfg.DCM.A)
    flag = [];
    for i = 1:3
        flag = [flag find(size(cfg.DCM.A{i})~=[Nareas, Nareas])];
    end
    if ~isempty(flag)
        error('Invalid A marix size!!');
    end
    DCM_low.A = cfg.DCM.A;
end
if ischar(cfg.DCM.B)
    DCM_low.B{1} = zeros(Nareas,Nareas); % B{1} : Effect dependent modulatory
elseif iscell(cfg.DCM.B)
    flag = [];
    for i = 1:length(cfg.DCM.B)
        flag = [flag find(size(cfg.DCM.A{i})~=[Nareas, Nareas])];
    end
    if ~isempty(flag)
        error('Invalid B marix size!!');
    end
    DCM_low.B = cfg.DCM.B;
end
DCM_low.C              = sparse(size(DCM_low.Lpos,2),0);
DCM_low.xU.X = sparse(1 ,0);
% C matrix is not used in CSD.
% Invert
% Spatial model
clear spm_erp_L

DCM_low.M.dipfit.model = DCM_low.options.model;
DCM_low.M.dipfit.type  = DCM_low.options.spatial;
DCM_low.options.DATA   = 1;
DCM_low                = spm_dcm_erp_data(DCM_low);             % data
DCM_low                = spm_dcm_erp_dipfit_jhs(DCM_low, 1);    % spatial model
Ns                 = length(DCM_low.A{1});                  % number of sources
% prior moments on parameters
[pE,pC]  = spm_dcm_neural_priors(DCM_low.A,DCM_low.B,DCM_low.C,DCM_low.options.model);
% augment with priors on spatial model
[pE,pC]  = spm_L_priors(DCM_low.M.dipfit,pE,pC); 
% augment with priors on endogenous inputs (neuronal) and noise
[pE,pC]  = spm_ssr_priors(pE,pC);
% initial states and equations of motion
[x,f]    = spm_dcm_x_neural(pE,DCM_low.options.model);

% create DCM
%--------------------------------------------------------------------------
DCM_low.M.IS = 'spm_csd_mtf';
DCM_low.M.FS = 'spm_fs_csd';
DCM_low.M.g  = 'spm_gx_erp';
DCM_low.M.f  = f;
DCM_low.M.x  = x;
DCM_low.M.n  = length(spm_vec(x));
DCM_low.M.pE = pE;
DCM_low.M.pC = pC;
DCM_low.M.hE = 6;
DCM_low.M.hC = 1/128; %1/64;
DCM_low.M.m  = Ns;

% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
DCM_low.M.u  = sparse(Ns,1);

%-Feature selection using principal components (U) of lead-field
%==========================================================================
 
% Spatial modes
%--------------------------------------------------------------------------
try
    DCM_low.M.U  = spm_dcm_eeg_channelmodes(DCM_low.M.dipfit,Nm);
end
 
% get data-features (in reduced eigenspace)
%==========================================================================
if DCM_low.options.DATA
    DCM_low  = spm_dcm_csd_data(DCM_low);
end
 
% scale data features (to a variance of about 8)
%--------------------------------------------------------------------------
ccf      = spm_csd2ccf(DCM_low.xY.y,DCM_low.xY.Hz);
scale    = max(spm_vec(ccf));
DCM_low.xY.y = spm_unvec(8*spm_vec(DCM_low.xY.y)/scale,DCM_low.xY.y);


% complete model specification and invert
%==========================================================================
Nm       = size(DCM_low.M.U,2);                    % number of spatial modes
DCM_low.M.l  = Nm;
DCM_low.M.Hz = DCM_low.xY.Hz;
DCM_low.M.dt = DCM_low.xY.dt;
 
% precision of noise: AR(1/2)
%--------------------------------------------------------------------------
y     = spm_fs_csd(DCM_low.xY.y,DCM_low.M);
for i = 1:length(y)
    n      = size(y{i},1);
    m      = size(y{i},2)*size(y{i},3);
    q      = spm_Q(1/2,n,1);
    Q{i,i} = kron(speye(m,m),q);
end
DCM_low.xY.Q  = spm_cat(Q);
DCM_low.xY.X0 = sparse(size(Q,1),0);


% Variational Laplace: model inversion
%==========================================================================
% test
DCM_low.M.pC.A{1}=DCM_low.M.pC.A{1}*10.0;
DCM_low.M.pC.A{2}=DCM_low.M.pC.A{2}*10.0;
DCM_low.M.pC.A{3}=DCM_low.M.pC.A{3}*10.0;
[Qp,Cp,Eh,F] = spm_nlsi_GN(DCM_low.M,DCM_low.xU,DCM_low.xY);


% Data ID
%--------------------------------------------------------------------------
try
    try
        ID = spm_data_id(feval(DCM_low.M.FS,DCM_low.xY.y,DCM_low.M));
    catch
        ID = spm_data_id(feval(DCM_low.M.FS,DCM_low.xY.y));
    end
catch
    ID = spm_data_id(DCM_low.xY.y);
end
 
 
% Bayesian inference {threshold = prior} NB Prior on A,B and C = exp(0) = 1
%==========================================================================
warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning('on', 'SPM:negativeVariance');
 
 
% predictions (csd) and error (sensor space)
%--------------------------------------------------------------------------
Hc  = spm_csd_mtf(Qp,DCM_low.M,DCM_low.xU);                      % prediction
Ec  = spm_unvec(spm_vec(DCM_low.xY.y) - spm_vec(Hc),Hc);     % prediction error
 
 
% predictions (source space - cf, a LFP from virtual electrode)
%--------------------------------------------------------------------------
M             = rmfield(DCM_low.M,'U'); 
M.dipfit.type = 'LFP';

M.U         = 1; 
M.l         = Ns;
qp          = Qp;
qp.L        = ones(1,Ns);             % set virtual electrode gain to unity
qp.b        = qp.b - 32;              % and suppress non-specific and
qp.c        = qp.c - 32;              % specific channel noise

[Hs Hz dtf] = spm_csd_mtf(qp,M,DCM_low.xU);
[ccf pst]   = spm_csd2ccf(Hs,DCM_low.M.Hz);
[coh fsd]   = spm_csd2coh(Hs,DCM_low.M.Hz);
DCM_low.dtf     = dtf;
DCM_low.ccf     = ccf;
DCM_low.coh     = coh;
DCM_low.fsd     = fsd;
DCM_low.pst     = pst;
DCM_low.Hz      = Hz;

 
% store estimates in DCM
%--------------------------------------------------------------------------
DCM_low.Ep = Qp;                   % conditional expectation
DCM_low.Cp = Cp;                   % conditional covariance
DCM_low.Pp = Pp;                   % conditional probability
DCM_low.Hc = Hc;                   % conditional responses (y), channel space
DCM_low.Rc = Ec;                   % conditional residuals (y), channel space
DCM_low.Hs = Hs;                   % conditional responses (y), source space
DCM_low.Ce = exp(-Eh);             % ReML error covariance
DCM_low.F  = F;                    % Laplace log evidence
DCM_low.ID = ID;                   % data ID
 
% and save
%--------------------------------------------------------------------------
DCM_low.options.Nmodes = Nm;


%% DCM High Alpha
DCM_high = [];
DCM_high.xY.Dfile = D.fullfile;
DCM_high.xY.modality = 'MEG';
DCM_high.name = ['DCM_' D.fname '_HighAlpha'];
DCM_high.val =D.val;
DCM_high.options.analysis = cfg.DCM.analysis;
DCM_high.options.model = cfg.DCM.model;
DCM_high.options.spatial = cfg.DCM.spatial;
DCM_high.options.trials = cfg.DCM.trials{2};
DCM_high.options.Tdcm = cfg.DCM.TimeWin;
DCM_high.options.Fdcm = cfg.DCM.FreqWin;
DCM_high.options.Nmodes = cfg.DCM.chanMode;
DCM_high.options.h = 1;
DCM_high.options.D = 1;
DCM_high.options.lock = 1;
DCM_high.options.multiC = 0;
DCM_high.options.symmetry = 0;
DCM_high.options.Rft = 5;
DCM_high = spm_dcm_erp_data(DCM_high, 0); % 0 for not to make erp
%spm_dcm_erp_results(DCM, 'Data');

% Prepare dipolefit. ( Set DCM.M . Only need for 'ECD' and 'IMG'. No need for 'LFP)
DCM_high.Lpos   = cfg.DCM.ROI;
DCM_high.Sname  = cfg.DCM.ROIName;
Nareas = size(DCM_high.Lpos,2);
% Spatial Model
DCM_high = spm_dcm_erp_dipfit_jhs(DCM_high);

% Prepare Connectivity Matrices ( Set A, B, C - No need C for 'CSD' )
if ischar(cfg.DCM.effect)
    DCM_high.xU.name = {}; % If there is some between trial effect, then you should specify xU.
    DCM_high.xU.X    = zeros(length(D.condlist),0); % Design Matrix
else
    if size(cfg.DCM.effect,1) == length(D.condlist)
        DCM_high.xU.name = cfg.DCM.effectName;
        DCM_high.xU.X    = cfg.DCM.effect;
    else
        warning('Size of Design Matrix does not match with the number of conditions');
        warning('By Default, we do not use Design matrix.');
        DCM_high.xU.name = {}; % If there is some between trial effect, then you should specify xU.
        DCM_high.xU.X    = zeros(length(D.condlist),0); % Design Matrix
    end
end
% Specify connectivity model ( Set as full model )
if ischar(cfg.DCM.A)
    DCM_high.A{1} = ones(Nareas, Nareas); % A{1} : Forward
    DCM_high.A{2} = ones(Nareas,Nareas); % A{2} : Backward
    DCM_high.A{3} = zeros(Nareas,Nareas); % A{3} : Modulatory
elseif iscell(cfg.DCM.A)
    flag = [];
    for i = 1:3
        flag = [flag find(size(cfg.DCM.A{i})~=[Nareas, Nareas])];
    end
    if ~isempty(flag)
        error('Invalid A marix size!!');
    end
    DCM_high.A = cfg.DCM.A;
end
if ischar(cfg.DCM.B)
    DCM_high.B{1} = zeros(Nareas,Nareas); % B{1} : Effect dependent modulatory
elseif iscell(cfg.DCM.B)
    flag = [];
    for i = 1:length(cfg.DCM.B)
        flag = [flag find(size(cfg.DCM.A{i})~=[Nareas, Nareas])];
    end
    if ~isempty(flag)
        error('Invalid B marix size!!');
    end
    DCM_high.B = cfg.DCM.B;
end
DCM_high.C              = sparse(size(DCM_high.Lpos,2),0);
DCM_high.xU.X = sparse(1 ,0);
% C matrix is not used in CSD.
% Invert
% Spatial model
clear spm_erp_L

DCM_high.M.dipfit.model = DCM_high.options.model;
DCM_high.M.dipfit.type  = DCM_high.options.spatial;
DCM_high.options.DATA   = 1;
DCM_high                = spm_dcm_erp_data(DCM_high);             % data
DCM_high                = spm_dcm_erp_dipfit_jhs(DCM_high, 1);    % spatial model
Ns                 = length(DCM_high.A{1});                  % number of sources
% prior moments on parameters
[pE,pC]  = spm_dcm_neural_priors(DCM_high.A,DCM_high.B,DCM_high.C,DCM_high.options.model);
% augment with priors on spatial model
[pE,pC]  = spm_L_priors(DCM_high.M.dipfit,pE,pC); 
% augment with priors on endogenous inputs (neuronal) and noise
[pE,pC]  = spm_ssr_priors(pE,pC);
% initial states and equations of motion
[x,f]    = spm_dcm_x_neural(pE,DCM_high.options.model);

% create DCM
%--------------------------------------------------------------------------
DCM_high.M.IS = 'spm_csd_mtf';
DCM_high.M.FS = 'spm_fs_csd';
DCM_high.M.g  = 'spm_gx_erp';
DCM_high.M.f  = f;
DCM_high.M.x  = x;
DCM_high.M.n  = length(spm_vec(x));
DCM_high.M.pE = pE;
DCM_high.M.pC = pC;
DCM_high.M.hE = 6;
DCM_high.M.hC = 1/256; %1/64;
DCM_high.M.m  = Ns;

% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
DCM_high.M.u  = sparse(Ns,1);

%-Feature selection using principal components (U) of lead-field
%==========================================================================
 
% Spatial modes
%--------------------------------------------------------------------------
try
    DCM_high.M.U  = spm_dcm_eeg_channelmodes(DCM_high.M.dipfit,Nm);
end
 
% get data-features (in reduced eigenspace)
%==========================================================================
if DCM_high.options.DATA
    DCM_high  = spm_dcm_csd_data(DCM_high);
end
 
% scale data features (to a variance of about 8)
%--------------------------------------------------------------------------
ccf      = spm_csd2ccf(DCM_high.xY.y,DCM_high.xY.Hz);
scale    = max(spm_vec(ccf));
DCM_high.xY.y = spm_unvec(8*spm_vec(DCM_high.xY.y)/scale,DCM_high.xY.y);


% complete model specification and invert
%==========================================================================
Nm       = size(DCM_high.M.U,2);                    % number of spatial modes
DCM_high.M.l  = Nm;
DCM_high.M.Hz = DCM_high.xY.Hz;
DCM_high.M.dt = DCM_high.xY.dt;
 
% precision of noise: AR(1/2)
%--------------------------------------------------------------------------
y     = spm_fs_csd(DCM_high.xY.y,DCM_high.M);
for i = 1:length(y)
    n      = size(y{i},1);
    m      = size(y{i},2)*size(y{i},3);
    q      = spm_Q(1/2,n,1);
    Q{i,i} = kron(speye(m,m),q);
end
DCM_high.xY.Q  = spm_cat(Q);
DCM_high.xY.X0 = sparse(size(Q,1),0);


% Variational Laplace: model inversion
%==========================================================================
% test
DCM_high.M.pC.A{1}=DCM_high.M.pC.A{1}*10.0;
DCM_high.M.pC.A{2}=DCM_high.M.pC.A{2}*10.0;
DCM_high.M.pC.A{3}=DCM_high.M.pC.A{3}*10.0;
[Qp,Cp,Eh,F] = spm_nlsi_GN(DCM_high.M,DCM_high.xU,DCM_high.xY);


% Data ID
%--------------------------------------------------------------------------
try
    try
        ID = spm_data_id(feval(DCM_high.M.FS,DCM_high.xY.y,DCM_high.M));
    catch
        ID = spm_data_id(feval(DCM_high.M.FS,DCM_high.xY.y));
    end
catch
    ID = spm_data_id(DCM_high.xY.y);
end
 
 
% Bayesian inference {threshold = prior} NB Prior on A,B and C = exp(0) = 1
%==========================================================================
warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning('on', 'SPM:negativeVariance');
 
 
% predictions (csd) and error (sensor space)
%--------------------------------------------------------------------------
Hc  = spm_csd_mtf(Qp,DCM_high.M,DCM_high.xU);                      % prediction
Ec  = spm_unvec(spm_vec(DCM_high.xY.y) - spm_vec(Hc),Hc);     % prediction error
 
 
% predictions (source space - cf, a LFP from virtual electrode)
%--------------------------------------------------------------------------
M             = rmfield(DCM_high.M,'U'); 
M.dipfit.type = 'LFP';

M.U         = 1; 
M.l         = Ns;
qp          = Qp;
qp.L        = ones(1,Ns);             % set virtual electrode gain to unity
qp.b        = qp.b - 32;              % and suppress non-specific and
qp.c        = qp.c - 32;              % specific channel noise

[Hs Hz dtf] = spm_csd_mtf(qp,M,DCM_high.xU);
[ccf pst]   = spm_csd2ccf(Hs,DCM_high.M.Hz);
[coh fsd]   = spm_csd2coh(Hs,DCM_high.M.Hz);
DCM_high.dtf     = dtf;
DCM_high.ccf     = ccf;
DCM_high.coh     = coh;
DCM_high.fsd     = fsd;
DCM_high.pst     = pst;
DCM_high.Hz      = Hz;

 
% store estimates in DCM
%--------------------------------------------------------------------------
DCM_high.Ep = Qp;                   % conditional expectation
DCM_high.Cp = Cp;                   % conditional covariance
DCM_high.Pp = Pp;                   % conditional probability
DCM_high.Hc = Hc;                   % conditional responses (y), channel space
DCM_high.Rc = Ec;                   % conditional residuals (y), channel space
DCM_high.Hs = Hs;                   % conditional responses (y), source space
DCM_high.Ce = exp(-Eh);             % ReML error covariance
DCM_high.F  = F;                    % Laplace log evidence
DCM_high.ID = ID;                   % data ID
 
% and save
%--------------------------------------------------------------------------
DCM_high.options.Nmodes = Nm;

