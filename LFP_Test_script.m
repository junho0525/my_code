%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LFP analysis code                                                       %
%    This LFP data is obtained by Jiye Choi and Prof. Yong Jeong in KAIST %
%    There is 2 mouse conditions and 3 age conditions (total 6 conditions)%
%        Tg : Transgenic mouse of Alzheimer's Disease                     %
%        Wt : Wild type mouse                                             %
%        Tg3, Tg7, Tg11, Wt3, Wt7, Wt11 : 3, 7 and 11 months aged mouse   %
%    The sample data is resting state data recorded from 8 ROIs(8channels)%
%        channel 1 : right anterior cingulate cortex (rACC)               %
%        channel 2 : left anterior cingulate cortex (lACC)                %
%        channel 3 : right anterior insula (rAI)                          %
%        channel 4 : left anterior insula (lAI)                           %
%        channel 5 : right retrosplenial cortex (rRSC)                    %
%        channel 6 : left retrosplenial cortex (lRSC)                     %
%        channel 7 : right basolateral amygdala (rBLA)                    %
%        channel 8 : left basolateral amygdala (lBLA)                     %
%    Sampling Frequency : 1600 Hz                                         %
%    Data Length : 10 seconds                                             %
%    This data is raw data, need to preprocess(60Hz notch filter, Detrend)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.09.18 11:47 - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Raw data
addpath('E:\my_code');
cd('11mo/WT/#4');
files = {'1-1.txt','2-1.txt','3-1.txt',...
    '4-1.txt','5-1.txt','6-1.txt',...
    '7-1.txt','8-1.txt','9-1.txt',...
    '10-1.txt'};

for session = 1:length(files)
    fid = fopen(files{session})
    Raw_d{session} = textscan(fid, '%f %f  %f  %f %f  %f %f  %f',20000);
end
fclose('all');
%% Check raw data
isdataok = 1;
for session = 1:length(Raw_d)
    if length(Raw_d{session})~=8
        isdataok = 0;
        error('This session has less than 8 channels')
    end
    for channel = 1:8
        if length(Raw_d{session}{channel})~=16000
            isdataok = 0;
            error ('This trial has more than 16000 time point');
        end
    end
end
if isdataok
    fprintf('Raw data is successfully loaded.');
end
%% Create fieldtrip type data(Concatenate all trials) 
for session = 1:10
    start(session) = ((session-1)*16000+1);
    epoch = [start(session): (start(session)+16000-1)];
    for chan = 1:8
        data(epoch,chan) = Raw_d{session}{chan};
    end
end

ftdata = [];
ftdata.trial{1} = data(:,1:8)';
ftdata.time{1} = (0:(160000-1))./1600;

% Channel define
% 1. right ACC  2.left ACC   3. right AI(Anterior Insula)   4. left AI
% 5. right RSC(Retrosplenial cortex)  6. left RSC 
% 7. right BLA(Basolateral Amygdala   8. left BLA 
chanlabel ={'rACC','lACC','rAI', 'lAI', 'rRSC', 'lRSC', 'rBLA','lBLA'};
ftdata.fsample = 1600;
ftdata.label = chanlabel(1:8);
ftdata.label = ftdata.label(:);
%% Preprocess(Downsample -> High pass filter-> Notch filter-> Low pass filter)
% Plot rawdata
figure;
subplot(4,2,1);
plot(ftdata.time{1}, ftdata.trial{1});
title('raw data timeseries(concatenated all trials)');
subplot(4,2,2);
plot((0:(160000-1))*(1600/160000),abs(fft(ftdata.trial{1}(1,:))));
title('raw data powerspectrum');

% Downample to 800Hz
cfg = [];
cfg.resamplefs = 800;
cfg.detrend = 'yes';
cfg.demean = 'yes';
datadown = ft_resampledata(cfg, ftdata);

subplot(4,2,3);
plot(datadown.time{1}, datadown.trial{1});
title('downsampled data timeseries (800 Hz)');
subplot(4,2,4);
plot((0:(80000-1))*(800/80000),abs(fft(datadown.trial{1}(1,:))));
title('downsampled data powerspectrum');

% Band pass filter (1, 200) Hz
cfg = [];
cfg.bpfilter ='yes';
cfg.bpfilttype = 'fir';
cfg.bpfreq = [0.5, 200]; %% must satisfy this condition ( freqmax*2 < fsample )
databpfir = ft_preprocessing(cfg, datadown);

subplot(4,2,5);
plot(databpfir.time{1}, databpfir.trial{1});
title('Band pass filtered data timeseries (1-200 Hz)');
subplot(4,2,6);
plot((0:(80000-1))*(800/80000),abs(fft(databpfir.trial{1}(1,:))));
title('Band pass filtered data powerspectrum');

% Notch filter(60Hz, 120Hz, 180Hz)
cfg = [];
cfg.dftfreq   = [60-1:(1/100):60+1 120-1:(1/100):120+1 ...
    180-1:(1/100):180+1]; % filter out 60 hz line noise
cfg.dftfilter = 'yes';
dataclean     = ft_preprocessing(cfg,databp);

subplot(4,2,7);
plot(dataclean.time{1}, dataclean.trial{1});
title('Notch filtered data timeseries (60, 120, 180 Hz)');
subplot(4,2,8);
plot((0:(80000-1))*(800/80000),abs(fft(dataclean.trial{1}(1,:))));
title('Notch pass filtered data powerspectrum');
%% Clean unnecessary data
clear data cfg chanlabel databp datadown epoch fid ftdata session
ftdata_concat = dataclean;
% separate each trials
ftdata = dataclean;
ftdata = rmfield(ftdata, 'sampleinfo');
ftdata.trial ={}; ftdata.time = {};
fsample = 800;
totsample = fsample * 10;
time = (0:(totsample-1))./fsample;
start = ((1:length(files))-1)*totsample +1;

for trl = 1:length(files)
    dataind = (0:(totsample-1))+start(trl);
    ftdata.trial{trl} = dataclean.trial{1}(:,dataind);
    ftdata.time{trl} = time;
end

clear dataclean
%% Create spm MEEG datafiles
% define the output file name
fname = 'LFP_WT_subj4';
fname_concat = 'LFP_WT_subj4_concat';

% Convert the ftdata struct to SPM M\EEG dataset
%--------------------------------------------------------------------------
D = spm_eeg_ft2spm(ftdata, fname);
D_concat = spm_eeg_ft2spm(ftdata_concat, fname_concat);

% Examples of providing additional information in a script
% ':' comes instead of an index vector and means that the command
% applies to all channels/all trials.
%--------------------------------------------------------------------------
D = type(D, 'single');             % Sets the dataset type
D = chantype(D, ':', 'LFP');        % Sets the channel type 
D = conditions(D, ':', 'WT-11');  % Sets the condition label

D_concat = type(D_concat, 'continuous');             % Sets the dataset type
D_concat = chantype(D_concat, ':', 'LFP');        % Sets the channel type 
D_concat = conditions(D_concat, ':', 'WT-11');  % Sets the condition label
% save
%--------------------------------------------------------------------------
save(D);
save(D_concat);

%% Run spDCM
if ~isdir('DCMs')
    mkdir('DCMs');
end
DCM = [];
DCM.xY.Dfile = D.fullfile;
DCM.xY.modality = 'LFP';
DCM.name = ['DCM_' D.fname];
%DCM.val =D.val;
DCM.options.analysis = 'CSD';
DCM.options.model = 'CMC';
DCM.options.spatial = 'LFP';
DCM.options.trials = 1;
DCM.options.Tdcm = [1 10000];
DCM.options.Fdcm = [1 60];
DCM.options.Nmodes = 8;
DCM.options.h = 1;
DCM.options.D = 1;
DCM.options.lock = 1;
DCM.options.multiC = 0;
DCM.options.symmetry = 0;
DCM.options.Rft = 5;
DCM = spm_dcm_erp_data(DCM, 0); % 0 for not to make erp
DCM.Lpos   = [];
DCM.Sname  = ftdata.label;
DCM.xU.X = [];
DCM.xU.name = {};
DCM.B ={};
DCM.A{1} = ones(8);
DCM.A{2} = ones(8);
DCM.A{3} = ones(8);
DCM = spm_dcm_csd(DCM);
%% Prepare dipolefit. ( Set DCM.M . Only need for 'ECD' and 'IMG'. No need for 'LFP)

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
