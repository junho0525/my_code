%% HCP MEG Data를 DCM으로 돌리기 위한 routine

%% HCP MEG Data는 이미 어느정도 Preprocessing이 되어있다.
%  HCP MEG Processing Pipeline을 통해, datacheck, baddata, icaclass를 통한 data
%  reduction이 이루어져 있는 상태이다.
%  Resting state analysis의 경우, rmegpreproc을, 사용하자.


%%  Load HCP MEG DATA 
%   나중에 배치로 돌리기위해 수정하자.
clear;
datasetPath = '/media/kuro/DATA1/HCPMEG';
workingPath = '/home/kuro/Codes/M_code/HCP_meeg';
%datasetPath = '/home/kuro/Codes/M_code/SPM_Examples/EEG_Single/';
cd (datasetPath);
apath = '100307_MEG/100307_MEG_Restin_preproc/100307/MEG/Restin/rmegpreproc';
filename = '100307_MEG_3-Restin_rmegpreproc.mat';
%filename = 'subject1.bdf';

fieldtripData = load(fullfile(datasetPath,apath,filename)); % MEG Data는 field trip data로 저장되어있다.
%fieldtripData = spm_eeg_convert([datasetPath,filename]);

if isfield( fieldtripData,'time')
    data = spm_eeg_ft2spm(fieldtripData, fullfile(datasetPath, strcat('spm_',filename)));
else
    data = spm_eeg_ft2spm(fieldtripData.data, fullfile(datasetPath, strcat('spm_',filename)));
end
%   [f,p]        = uigetfile('*.mat','please select data file');
%   file이 안열릴때는 위의 코드를 이용해 유저가 직접 파일을 열수있도록 하자.
%%  Resting State Data Vectorize
%   왠지는 모르겠으나, HCP데이터는 하나의 time-series data를 여러 trial로 쪼갠다(resting임에도.)
%   따라서, 이를 하나의 time-series로서 재구성하는 과정.


%%  SPM MEG Preprocessing
%   SPM EEG 예제의 Preprocessing 절차를 그대로 따라간다.
%   High-pass filter -> Downsampling -> Low-pass filter
%    (Cutoff = 0.1)   (down to 200Hz)    (Cutoff = 30)
%   다만 Slow wave의 경우, 40Hz이상을 보기도 하기때문에, 이경우 Low-pass filter의 
%   Cutoff를 더 높게 조정해주면 된다.
%   filter의 경우, spm_eeg_filter를 사용한다.
%       FORMAT D = spm_eeg_filter(S)
%
%       S           - input structure
%        Fields of S:
%         S.D       - MEEG object or filename of M/EEG mat-file
%
%         S.band    - filterband [low|high|bandpass|stop]
%         S.freq    - cutoff frequency(-ies) [Hz]
%
%        Optional fields:
%         S.type    - filter type [default: 'butterworth']
%                         'butterworth': Butterworth IIR filter
%                         'fir':         FIR filter (using MATLAB fir1 function)
%         S.order   - filter order [default: 5 for Butterworth]
%         S.dir     - filter direction [default: 'twopass']
%                         'onepass'         forward filter only
%                         'onepass-reverse' reverse filter only, i.e. backward in time
%                         'twopass'         zero-phase forward and reverse filter
%         S.prefix  - prefix for the output file [default: 'f']
%
%        D           - MEEG object (also written to disk)

%%
%   High-pass Filter
S.D = data;
S.band = 'high';
S.freq = 0.1;
S.prefix = 'f_';
data = spm_eeg_filter(S);
clear S;

%%
%   Downsampling
S.D = data;
S.method = 'fft'; % 'resample' (default), 'decimate', 'downsample', 'fft'
S.fsample_new = 200;
data = spm_eeg_downsample(S);
clear S;

%%
%   Low-pass Filter (Result File will have 'fdf_spm_' prefix
S.D = data;
S.band = 'low';
S.freq = 30;
data = spm_eeg_filter(S);
clear S;

%%
%   Epoching Step - Trial들간의 gap을 메꿔주어 single trial처럼 만들어준다.
%   현재는 type이 single로 표시되어 있기때문에 안돌아간다. HCP에서 datacheck, baddata, icaclass가
%   정확히 어떠한 preprocessing 과정인지 확실히 알아야 할필요가 있을듯 하다.
% S.D = data;

%%
%   Averaging
S.D = data;
data = spm_eeg_average(S);
clear S;

%%
%   MEG Signal Visualization
spm_eeg_review(data);


%%
%--------------------------------------------------------------------------
%   DCM File 만들기
%sp = '/home/kuro/Codes/M_code/VSDI_DCM_by_KSJung/';
%cd(sp)
%DCMfile = fullfile(sp,'Before_invert','DCM_dLFP_before_invert.mat');
%load(DCMfile,'DCM') % DCM파일을 불러와본다.
name = fullfile(workingPath,strcat('DCM_',filename));
DCM.xY.Dfile = fullfile(datasetPath, strcat('mfdf_spm_',filename));

%%
clear spm_erp_L
dt = 5; % ms; sampling rate = 200 Hz
DCM.options.analysis  = 'CSD';
DCM.options.spatial   = 'IMG';
DCM.options.model     = 'CMC';
DCM.options.Tdcm      = [1,2000];     % [start end] time window in ms
DCM.options.Fdcm      = [1,200];     % [start end] Frequency window in Hz
DCM.options.Nmodes    = 8;      % # of feature
DCM.options.DATA      = 1;      % y/n preprocessed data
DCM.options.trials    = [1];    % # of trials
%DCM.options.onset     = 70;     % onset time of the first stimulus (in ms)
%DCM.options.dur       = 1;      % duration
DCM.options.D         = 1;      % Down-sampling interval

%DCM = rmfield(DCM,'xU');
% DCM.xU.X              = [1];
% DCM.xU.name           = {'Stim'};
% DCM.xU.u              = 
% DCM.xU.dt             = 


%% Specify DCM
% Filename and options
%--------------------------------------------------------------------------
try, DCM.name;                      catch, DCM.name = name;      end
try, model   = DCM.options.model;   catch, model    = 'NMM';     end
try, spatial = DCM.options.spatial; catch, spatial  = 'LFP';     end
try, Nm      = DCM.options.Nmodes;  catch, Nm       = 8;         end
try, DATA    = DCM.options.DATA;    catch, DATA     = 1;         end
DCM.options.Nmodes = Nm;
DCM.name = name;
%% Location Priors for dipoles
DCM.Lpos  = [[-42; -22; 7] [46; -14; 8] [-61; -32; 8] [59; -25; 8] [46; 20; 8]];
DCM.Sname = {'left AI', 'right A1', 'left STG', 'right STG', 'right IFG'};
Nareas    = size(DCM.Lpos,2);

%% DCM.M.dipfit 설정
DCM.M.dipfit.location = 0; % not allow source location change
DCM.M.dipfit.Lpos = DCM.Lpos;
DCM.M.dipfit.Nm = Nm;
DCM.M.dipfit.Ns = Nareas;
DCM.M.dipfit.Nc = size(data.chanlabels,2);
DCM.M.dipfit.model = model;
DCM.M.dipfit.type  = spatial;
DCM.M.dipfit.silent_source = {}; % silent source for CSD 'IMG'

%%
% Spatial model
%==========================================================================

DCM  = spm_dcm_erp_data(DCM);                   
DCM  = spm_dcm_erp_dipfit(DCM, 1);

if DATA
    DCM = spm_dcm_csd_data(DCM);
    DCM = spm_dcm_csd_dipfit(DCM, 1);
end


if DATA
    DCM  = spm_dcm_erp_data(DCM);                   % data
    DCM  = spm_dcm_erp_dipfit(DCM, 1);              % spatial model
end
Ns   = length(DCM.A{1});                            % number of sources
DCM.M.dipfit.Ns    = Ns; % KSJ





%%
%DCM = rmfield(DCM,{'A','B', 'Sname'});
DCM.A{1,1}            = ones(4,4);
DCM.B{1,1}            = ones(4,4);
DCM.C                 = zeros(4,1);
DCM.Sname{1,1}        = 'Full model'; %source name
DCM                   = rmfield(DCM,'xY');
DCM.xY.modality       = 'LFP';
DCM.xY.name           = {'Full model'};
DCM.xY.Time           = [0:dt:2000]; % all scans

%% Generate MEEG data from an abitrary data
%--------------------------------------------------------------------------
% Examples of providing additional information in a script
% [] comes instead of an index vector and means that the command
% applies to all channels/all trials.
%--------------------------------------------------------------------------
% Select Channels
%
%data_selected.trial{1} = data ([15, 100, 123,237], :,1);
%data_selected.time{1} = [0:dt:2000];

%data_selected.fsample = 200;
%data_selected.label = fieldtripData.data.label([15,100, 123, 237]);

%D = spm_eeg_ft2spm(data_selected, 'MEEG_converted_data');

%D = chantype(D, ':', 'LFP');                   % Sets the channel type 
%D = conditions(D, 1, 'Condition 1');  % Sets the condition label



% save 
%--------------------------------------------------------------------------
%save(D);
%clear data_selected;


%%
% Design model and exogenous inputs
%==========================================================================
if ~isfield(DCM,'xU'),   DCM.xU.X = sparse(1 ,0); end
if ~isfield(DCM.xU,'X'), DCM.xU.X = sparse(1 ,0); end
if ~isfield(DCM,'C'),    DCM.C    = sparse(Ns,0); end
if isempty(DCM.xU.X),    DCM.xU.X = sparse(1 ,0); end
if isempty(DCM.xU.X),    DCM.C    = sparse(Ns,0); end

% Neural mass model
%==========================================================================
 
% prior moments on parameters
%--------------------------------------------------------------------------
[pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,model);
  
% check to see if neuronal priors have already been specified
%--------------------------------------------------------------------------
try
    if spm_length(DCM.M.pE) == spm_length(pE);
        pE = DCM.M.pE;
        pC = DCM.M.pC;
        fprintf('Using existing priors\n')
    end
end
%%
% augment with priors on spatial model
%--------------------------------------------------------------------------
[pE,pC]  = spm_L_priors(DCM.M.dipfit,pE,pC);

%% 
% augment with priors on endogenous inputs (neuronal) and noise
%--------------------------------------------------------------------------
[pE,pC]  = spm_ssr_priors(pE,pC);

try
    if spm_length(DCM.M.pE) == spm_length(pE);
        pE = DCM.M.pE;
        pC = DCM.M.pC;
        fprintf('Using existing priors\n')
    end
end

%%
% initial states and equations of motion
%--------------------------------------------------------------------------
[x,f]    = spm_dcm_x_neural(pE,model);

%%
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
DCM.M.hE = 6;             % Hyper Prior의 Expectation 값을 이렇게 주는데 이유는 모르겠다.
DCM.M.hC = 1/64;
DCM.M.m  = Ns;

% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
DCM.M.u  = sparse(Ns,1);

%%
%-Feature selection using principal components (U) of lead-field
%==========================================================================
 
% Spatial modes
%--------------------------------------------------------------------------
try
    DCM.M.U  = spm_dcm_eeg_channelmodes(DCM.M.dipfit,Nm);
end

%%
% get data-features (in reduced eigenspace)
%==========================================================================
if DATA
    DCM  = spm_dcm_csd_data(DCM);
end

% scale data features (to a variance of about 8)
%--------------------------------------------------------------------------
ccf      = spm_csd2ccf(DCM.xY.y,DCM.xY.Hz);
scale    = max(spm_vec(ccf));
DCM.xY.y = spm_unvec(8*spm_vec(DCM.xY.y)/scale,DCM.xY.y);

%%
% complete model specification and invert
%==========================================================================
Nm       = size(DCM.M.U,2);                    % number of spatial modes
DCM.M.l  = Nm;
DCM.M.Hz = DCM.xY.Hz;
DCM.M.dt = DCM.xY.dt;

%%
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

save(DCM.name, 'DCM', spm_get_defaults('mat.format'));
%%
% Variational Laplace: model inversion
%==========================================================================
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
 
%% 
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
 
save(DCM.name, 'DCM', spm_get_defaults('mat.format'));
%%
function show_meeg(meeg_data)
    if ~isa(meeg_data, 'meeg')
        return
    end
end