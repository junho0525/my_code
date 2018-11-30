function [D,DCM] = mnet_hcp_meeg_rest_spDCM_erpmodel(config,D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEEG prprocessing for HCP dataset by using SPM.                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.05.14 12:12 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check Input
if nargin == 1
    % These five fields must be specified
    if ~isfield(config, 'subjectID'),error('Invlid Subject ID!! There is no subjectID');end
    if (~isfield(config, 'subjectPath')||~isdir(config.subjectPath));error('Invalid Subject Path!! There is no such folder');end
    if (~isfield(config, 'savePath')||~isdir(config.savePath)),error('Invalid Save Path!! There is no such folder');end
    if (~isfield(config, 'sessionID')),error('Invalid Session ID!! There is no session ID field!!');end
    %% Set Path
    rMEGPath = fullfile(config.subjectPath,[config.subjectID '_MEG_Restin_preproc'],config.subjectID,'MEG','Restin','rmegpreproc');
    rMEGRawData=[config.subjectID '_MEG_' num2str(config.sessionID) '-Restin_rmegpreproc'];
    %% Load MEEG Datas
    try
        D = spm_eeg_load(fullfile(config.savePath,['affdspm_', rMEGRawData]));
        fieldtripData = load(fullfile(rMEGPath,rMEGRawData));
    catch
        error(['There is no preprocessed file!! filename : ' fullfile(config.savePath,['affdspm_', rMEGRawData])]);
    end
end
if ~(isfield(D,'inv') &&isfield(D.inv{1},'forward')&&isfield(D.inv{1},'datareg')&&isfield(D.inv{1},'mesh'))
    error('There is no forward model!!');
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
DCM.options.Fdcm = config.FreqWin; % [1 40]
DCM.options.Nmodes = config.chanMode;
DCM.options.h = 1;
DCM.options.D = 1;
DCM.options.lock = 1;
DCM.options.multiC = 0;
DCM.options.symmetry = 0;
DCM.options.Rft = 5;
DCM = spm_dcm_erp_data(DCM, 0); % 0 for not to make erp
%spm_dcm_erp_results(DCM, 'Data');

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
DCM.C              = sparse(size(DCM.Lpos,2),0);
DCM.xU.X = sparse(1 ,0);
% C matrix is not used in CSD.
%% Invert
% Spatial model
clear spm_erp_L

DCM.M.dipfit.model = DCM.options.model;
DCM.M.dipfit.type  = DCM.options.spatial;
DCM.options.DATA   = 1;
DCM                = spm_dcm_erp_data(DCM);             % data
DCM                = spm_dcm_erp_dipfit_jhs(DCM, 1);    % spatial model
Ns                 = length(DCM.A{1});                  % number of sources
% prior moments on parameters
[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
% augment with priors on spatial model
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC); 
% augment with priors on endogenous inputs (neuronal) and noise
[pE,pC]  = spm_ssr_priors(pE,pC);
% initial states and equations of motion
[x,f]    = spm_dcm_x_neural(pE,DCM.options.model);

% create DCM
%--------------------------------------------------------------------------
DCM.M.IS = 'spm_csd_mtf';
DCM.M.FS = 'spm_fs_csd';
DCM.M.g  = 'spm_gx_erp';
DCM.M.f  = f;
DCM.M.x  = x;
DCM.M.n  = length(spm_vec(x));
DCM.M.pE = pE;
DCM.M.pC = pC;
DCM.M.hE = 6;
DCM.M.hC = 1/128; %1/64;
DCM.M.m  = Ns;

% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
DCM.M.u  = sparse(Ns,1);

%-Feature selection using principal components (U) of lead-field
%==========================================================================
 
% Spatial modes
%--------------------------------------------------------------------------
try
    DCM.M.U  = spm_dcm_eeg_channelmodes(DCM.M.dipfit,Nm);
end
 
% get data-features (in reduced eigenspace)
%==========================================================================
if DCM.options.DATA
    DCM  = spm_dcm_csd_data(DCM);
end
 
% scale data features (to a variance of about 8)
%--------------------------------------------------------------------------
ccf      = spm_csd2ccf(DCM.xY.y,DCM.xY.Hz);
scale    = max(spm_vec(ccf));
DCM.xY.y = spm_unvec(8*spm_vec(DCM.xY.y)/scale,DCM.xY.y);


% complete model specification and invert
%==========================================================================
Nm       = size(DCM.M.U,2);                    % number of spatial modes
DCM.M.l  = Nm;
DCM.M.Hz = DCM.xY.Hz;
DCM.M.dt = DCM.xY.dt;
 
% precision of noise: AR(1/2)
%--------------------------------------------------------------------------
y     = spm_fs_csd(DCM.xY.y,DCM.M);
for i = 1:length(y)
    n      = size(y{i},1);
    m      = size(y{i},2)*size(y{i},3);
    q      = spm_Q(1/2,n,1);
    Q{i,i} = kron(speye(m,m),q);
end
DCM.xY.Q  = spm_cat(Q);
DCM.xY.X0 = sparse(size(Q,1),0);


% Variational Laplace: model inversion
%==========================================================================
% test
DCM.M.pC.A{1}=DCM.M.pC.A{1}*10.0;
DCM.M.pC.A{2}=DCM.M.pC.A{2}*10.0;
DCM.M.pC.A{3}=DCM.M.pC.A{3}*10.0;
[Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);


% Data ID
%--------------------------------------------------------------------------
try
    try
        ID = spm_data_id(feval(DCM.M.FS,DCM.xY.y,DCM.M));
    catch
        ID = spm_data_id(feval(DCM.M.FS,DCM.xY.y));
    end
catch
    ID = spm_data_id(DCM.xY.y);
end
 
 
% Bayesian inference {threshold = prior} NB Prior on A,B and C = exp(0) = 1
%==========================================================================
warning('off','SPM:negativeVariance');
dp  = spm_vec(Qp) - spm_vec(pE);
Pp  = spm_unvec(1 - spm_Ncdf(0,abs(dp),diag(Cp)),Qp);
warning('on', 'SPM:negativeVariance');
 
 
% predictions (csd) and error (sensor space)
%--------------------------------------------------------------------------
Hc  = spm_csd_mtf(Qp,DCM.M,DCM.xU);                      % prediction
Ec  = spm_unvec(spm_vec(DCM.xY.y) - spm_vec(Hc),Hc);     % prediction error
 
 
% predictions (source space - cf, a LFP from virtual electrode)
%--------------------------------------------------------------------------
M             = rmfield(DCM.M,'U'); 
M.dipfit.type = 'LFP';

M.U         = 1; 
M.l         = Ns;
qp          = Qp;
qp.L        = ones(1,Ns);             % set virtual electrode gain to unity
qp.b        = qp.b - 32;              % and suppress non-specific and
qp.c        = qp.c - 32;              % specific channel noise

[Hs Hz dtf] = spm_csd_mtf(qp,M,DCM.xU);
[ccf pst]   = spm_csd2ccf(Hs,DCM.M.Hz);
[coh fsd]   = spm_csd2coh(Hs,DCM.M.Hz);
DCM.dtf     = dtf;
DCM.ccf     = ccf;
DCM.coh     = coh;
DCM.fsd     = fsd;
DCM.pst     = pst;
DCM.Hz      = Hz;

 
% store estimates in DCM
%--------------------------------------------------------------------------
DCM.Ep = Qp;                   % conditional expectation
DCM.Cp = Cp;                   % conditional covariance
DCM.Pp = Pp;                   % conditional probability
DCM.Hc = Hc;                   % conditional responses (y), channel space
DCM.Rc = Ec;                   % conditional residuals (y), channel space
DCM.Hs = Hs;                   % conditional responses (y), source space
DCM.Ce = exp(-Eh);             % ReML error covariance
DCM.F  = F;                    % Laplace log evidence
DCM.ID = ID;                   % data ID
 
% and save
%--------------------------------------------------------------------------
DCM.options.Nmodes = Nm;
