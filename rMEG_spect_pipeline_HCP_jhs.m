%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Resting State MEG Frequency Analysis Pipeline Using SPM          %
%                                                                         %
% This pipeline code is for HCP dataset                                   %
%     Example data - 100307_MEG_3-Restin_rmegpreproc                      %
%                                                                         %
% Source localization with different frequency band.(eg. from delta to    %
% gamma)                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2018.01.26 xx:xx - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load HCP Data
%

clear;
close all;
try 
    D = spm_eeg_load('spm_100307_MEG_3-Restin_rmegpreproc');
    filename = '100307_MEG_3-Restin_rmegpreproc';
    fieldtripData = load(filename); % HCP MEG Data is field trip structure data
catch
    filename = '100307_MEG_3-Restin_rmegpreproc';
    fieldtripData = load(filename); % HCP MEG Data is field trip structure data
    if isfield( fieldtripData,'time') % Convert Field trip data to spm meeg object and pefixed with 'spm_'
        D = spm_eeg_ft2spm(fieldtripData, ['spm_', filename]);
    else
        D = spm_eeg_ft2spm(fieldtripData.data, ['spm_', filename]);
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SPM MEG Preprocessing
%
%   This step follows SPM EEG Example(Ch.40 of SPM manual) Preprocessing
%   steps
%   
%   Preprocessing steps are shown below
%
%   High-pass filter   ->   Downsampling      ->      Low-pass filter
%  (Cutoff = 0.1Hz)    (down to 200 samples/sec)     (Cutoff = 30Hz)
%
%   But for some cases, over 40Hz signals may have valuable information,
%   user should change low-pass filter cutoff value to higher than 40Hz.
%   
%   for Low/High - pass filter, we use spm_eeg_filter() function.
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaD = D; %delta wave - 2-4 Hz
thetaD = D; %theta wave - 5-7 Hz
alphaD = D; %alpha wave - 8-12 Hz
betaD = D; %beta wave - 15-29 Hz
gamma1D = D; %gamma 1 wave - 30-60 Hz
gamma2D = D; %gamma 2 wave - 60-90 Hz


%   High-pass Filter
S.D = D;
S.band = 'high';
S.freq = 2;
S.prefix = 'f_';
D = spm_eeg_filter(S);
clear S;
%   Downsampling
S.D = D;
S.method = 'fft'; % 'resample' (default), 'decimate', 'downsample', 'fft'
S.fsample_new = 200;
D = spm_eeg_downsample(S);
clear S;
%   Low-pass Filter (Result File will have 'fdf_spm_' prefix
S.D = D;
S.band = 'low';
S.freq = 90;
D = spm_eeg_filter(S);
clear S;

% spectrum data
%   Delta ( 2 - 4 Hz )
%   High-pass Filter
S.D = deltaD;
S.band = 'high';
S.freq = 2;
S.prefix = 'f_delta';
deltaD = spm_eeg_filter(S);
clear S;
%   Downsampling
S.D = deltaD;
S.method = 'fft'; % 'resample' (default), 'decimate', 'downsample', 'fft'
S.fsample_new = 200;
deltaD = spm_eeg_downsample(S);
clear S;
%   Low-pass Filter (Result File will have 'fdf_spm_' prefix
S.D = deltaD;
S.band = 'low';
S.freq = 4;
deltaD = spm_eeg_filter(S);
clear S;

%   Theta ( 5 - 7 Hz )
%   High-pass Filter
S.D = thetaD;
S.band = 'high';
S.freq = 5;
S.prefix = 'f_theta';
thetaD = spm_eeg_filter(S);
clear S;
%   Downsampling
S.D = thetaD;
S.method = 'fft'; % 'resample' (default), 'decimate', 'downsample', 'fft'
S.fsample_new = 200;
thetaD = spm_eeg_downsample(S);
clear S;
%   Low-pass Filter (Result File will have 'fdf_spm_' prefix
S.D = thetaD;
S.band = 'low';
S.freq = 7;
thetaD = spm_eeg_filter(S);
clear S;

%   Alpha ( 8 - 12 Hz )
%   High-pass Filter
S.D = alphaD;
S.band = 'high';
S.freq = 8;
S.prefix = 'f_alpha';
alphaD = spm_eeg_filter(S);
clear S;
%   Downsampling
S.D = alphaD;
S.method = 'fft'; % 'resample' (default), 'decimate', 'downsample', 'fft'
S.fsample_new = 200;
alphaD = spm_eeg_downsample(S);
clear S;
%   Low-pass Filter (Result File will have 'fdf_spm_' prefix
S.D = alphaD;
S.band = 'low';
S.freq = 12;
alphaD = spm_eeg_filter(S);
clear S;

%   Beta ( 15 - 29 Hz )
%   High-pass Filter
S.D = betaD;
S.band = 'high';
S.freq = 15;
S.prefix = 'f_beta';
betaD = spm_eeg_filter(S);
clear S;
%   Downsampling
S.D = betaD;
S.method = 'fft'; % 'resample' (default), 'decimate', 'downsample', 'fft'
S.fsample_new = 200;
betaD = spm_eeg_downsample(S);
clear S;
%   Low-pass Filter (Result File will have 'fdf_spm_' prefix
S.D = betaD;
S.band = 'low';
S.freq = 29;
betaD = spm_eeg_filter(S);
clear S;

%   Gamma1 ( 30 - 60 Hz )
%   High-pass Filter
S.D = gamma1D;
S.band = 'high';
S.freq = 30;
S.prefix = 'f_gamma1';
gamma1D = spm_eeg_filter(S);
clear S;
%   Downsampling
S.D = gamma1D;
S.method = 'fft'; % 'resample' (default), 'decimate', 'downsample', 'fft'
S.fsample_new = 200;
gamma1D = spm_eeg_downsample(S);
clear S;
%   Low-pass Filter (Result File will have 'fdf_spm_' prefix
S.D = gamma1D;
S.band = 'low';
S.freq = 60;
gamma1D = spm_eeg_filter(S);
clear S;

%   Gamma2 ( 60 - 90 Hz )
%   High-pass Filter
S.D = gamma2D;
S.band = 'high';
S.freq = 60;
S.prefix = 'f_gamma2';
gamma2D = spm_eeg_filter(S);
clear S;
%   Downsampling
S.D = gamma2D;
S.method = 'fft'; % 'resample' (default), 'decimate', 'downsample', 'fft'
S.fsample_new = 200;
gamma2D = spm_eeg_downsample(S);
clear S;
%   Low-pass Filter (Result File will have 'fdf_spm_' prefix
S.D = gamma2D;
S.band = 'low';
S.freq = 90;
gamma2D = spm_eeg_filter(S);
clear S;

%% Check bandpass filter
Fs = D.fsample;
L = size(D,2);
f = Fs*(0:(L/2))/L;

Y = fft(D(50,:,100));

P2 = abs(Y/L);
P1 = P2(1:(L/2+1));
P1(2:end-1) = 2*P1(2:end-1);

figure;
plot(f,P1);
hold on;

Y = fft(deltaD(50,:,100));
P2 = abs(Y/400);
P1 = P2(1:201);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1);

Y = fft(thetaD(50,:,100));
P2 = abs(Y/400);
P1 = P2(1:201);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1);

Y = fft(alphaD(50,:,100));
P2 = abs(Y/400);
P1 = P2(1:201);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1);

Y = fft(betaD(50,:,100));
P2 = abs(Y/400);
P1 = P2(1:201);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1);

Y = fft(gamma1D(50,:,100));
P2 = abs(Y/400);
P1 = P2(1:201);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1);

Y = fft(gamma2D(50,:,100));
P2 = abs(Y/400);
P1 = P2(1:201);
P1(2:end-1) = 2*P1(2:end-1);
plot(f,P1);
hold off;

%% Check Using Spectrogram
figure;
data = D(50, :, 100);
spectrogram(data, 32, 30, 256, Fs, 'yaxis');
figure;
data = alphaD(50, :, 100);
spectrogram(data, 32, 30, 256, Fs, 'yaxis'); 
%% Specify Forward Model
%
%
%
%%%opengl hardware;


return_path = pwd;

% ensure that the time and date of execution are not stored in the provenance information
%%%global ft_default
%%%ft_default.trackcallinfo = 'no';

hcp_read_matlab('100307_MEG_3-Restin_icaclass_vs.mat'); % This file is originally in Restin/icaclass directory of Restin_preproc

% read the source and volume conduction model from current dir with
% outputs of previous pipelines
fid = fopen('100307_MEG_anatomy_transform.txt', 'rt');
strFid=fread(fid,[1 inf], '*char');
eval(strFid);
fclose(fid);
clear fid;
clear strFid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify sourcemodel, headmodel, sensor and transform if it necessary
hcp_read_matlab('100307_MEG_anatomy_sourcemodel_2d.mat'); % This file is originally in anatomy directory
sourcemodel2d=ft_convert_units(sourcemodel2d, 'mm');
sourcemodelsubj = sourcemodel2d;

hcp_read_matlab('100307_MEG_anatomy_headmodel.mat');
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

cfg = [];
cfg.vol = headmodel;
cfg.grid = sourcemodelsubj;
cfg.grad = grad;
cfg.channel = channels;
cfg.normalize = 'yes';
cfg.reducerank = 2;

clear val tlck subjectid sourcemodelsubj resultprefix options mixing megpath;
clear rmeg mainpath i headmodel grad ft_default experimentid comp_class channels anatpath;

% specify forward model (Bring forward model information 
inv = [];
inv{1}.forward = [];
inv{1}.method='Imaging';
inv{1}.forward.toMNI = transform.bti2spm;
%
inv{1}.forward.fromMNI = transform.spm2bti;
inv{1}.forward.voltype='Single Shell';
inv{1}.forward.vol = ft_convert_units(cfg.vol,'m');
inv{1}.forward.modality = 'MEG';
inv{1}.forward.siunits=1;
inv{1}.forward.mesh=[];
inv{1}.forward.mesh_correction=[];
inv{1}.forward.mesh.face = cfg.grid.tri;
grid = ft_convert_units(cfg.grid,'m');
inv{1}.forward.mesh.vert = grid.pos;
inv{1}.forward.sensors = [];
inv{1}.forward.sensors = ft_convert_units(cfg.grad,'m');
inv{1}.mesh.tess_mni.face = inv{1}.forward.mesh.face;
grid = ft_convert_units(cfg.grid,'mm');
inv{1}.mesh.tess_mni.vert = spm_eeg_inv_transform_points(inv{1}.forward.toMNI,grid.pos); % tess_mni should have mni coordinate
Datareg(1).modality = 'MEG';
inv{1}.datareg(1) = Datareg(1);
inv{1}.forward.toMNI = [inv{1}.forward.toMNI(:, 1:3)*1000 inv{1}.forward.toMNI(:,4)];
inv{1}.forward.fromMNI = eye(4)/inv{1}.forward.toMNI;
D.inv = inv;
deltaD.inv = inv;
thetaD.inv = inv;
alphaD.inv = inv;
betaD.inv = inv;
gamma1D.inv = inv;
gamma2D.inv = inv;
clear Datareg grid inv;
% check meshes and co-registration
spm_eeg_inv_checkforward(D); 
% check meshes using fieldtrip
% figure;
% ft_plot_sens(D.inv{1}.forward.sensors);
% ft_plot_mesh(sourcemodel2d);
% ft_plot_vol(D.inv{1}.forward.vol);
%% Invert SPM Example
%  spm_cfg_eeg_inv_invert.m -> function run_inversion(job)
%  Specify job
%

job = [];
job.D = {D}; % job.D is cell type
job.val = 1; %inversion index
job.whatconditions.all = 1;
% Specify Inverting Methods.
% Use 'LOR', 'GS', 'IID', 'EEB' for Methods 'COH', 'MSP(GS)', 'IID', 'EBB'
job.isstandard.custom.invtype='GS'; 
job.isstandard.custom.woi = [0, 600];
job.isstandard.custom.foi = [0, 40]; % simulated data has 10, 20 Hz sources, so 0 to 40 Hz is enough to cover source activity
job.isstandard.custom.hanning = 1;
job.isstandard.custom.priors.priorsmask{1} = '';
job.isstandard.custom.priors.space = 1;
job.isstandard.custom.restrict.locs = [];
job.isstandard.custom.restrict.radius = 32;
job.isstandard.custom.restrict.mask{1}='';
job.modality{1} = 'All';


%
% Run Invert
% spm_cfg_eeg_inv_invert code
%

inverse = [];
if isfield(job.whatconditions, 'condlabel')
    inverse.trials = job.whatconditions.condlabel;
end

if isfield(job.isstandard, 'custom')
    inverse.type = job.isstandard.custom.invtype;
    inverse.woi  = fix([max(min(job.isstandard.custom.woi), 1000*D.time(1)) min(max(job.isstandard.custom.woi), 1000*D.time(end))]);
    inverse.Han  = job.isstandard.custom.hanning;
    inverse.lpf  =  fix(min(job.isstandard.custom.foi));
    inverse.hpf  =  fix(max(job.isstandard.custom.foi));
    
    P = char(job.isstandard.custom.priors.priorsmask);
    if ~isempty(P)        
        [p,f,e] = fileparts(P);
        switch lower(e)
            case '.gii'
                g = gifti(P);
                inverse.pQ = cell(1,size(g.cdata,2));
                for i=1:size(g.cdata,2)
                    inverse.pQ{i} = double(g.cdata(:,i));
                end
            case '.mat'
                load(P);
                inverse.pQ = pQ;
            case {'.img', '.nii'}
                S.D = D;
                S.fmri = P;
                S.space = job.isstandard.custom.priors.space;
                D = spm_eeg_inv_fmripriors(S);
                inverse.fmri = D.inv{D.val}.inverse.fmri;
                load(inverse.fmri.priors);
                inverse.pQ = pQ;
            otherwise
                error('Unknown file type.');
        end
    end
    
    if ~isempty(job.isstandard.custom.restrict.locs)
        inverse.xyz = job.isstandard.custom.restrict.locs;
        inverse.rad = job.isstandard.custom.restrict.radius;
    end
    
    P = char(job.isstandard.custom.restrict.mask);
    if ~isempty(P)
        inverse.mask = P;
    end
end

[mod, list] = modality(D, 1, 1);
if strcmp(job.modality{1}, 'All')
    inverse.modality  = list;
else
    inverse.modality  = intersect(list, job.modality);
end

if numel(inverse.modality) == 1
    inverse.modality = inverse.modality{1};
end
save(D);
D = {};

for i = 1:numel(job.D)
    D{i} = spm_eeg_load(job.D{i});
    
    D{i}.val = job.val;
    
    D{i}.con = 1;
    
    if ~isfield(D{i}, 'inv')
        error('Forward model is missing for subject %d.', i);
    elseif  numel(D{i}.inv)<D{i}.val || ~isfield(D{i}.inv{D{i}.val}, 'forward')
        if D{i}.val>1 && isfield(D{i}.inv{D{i}.val-1}, 'forward')
            D{i}.inv{D{i}.val} = D{i}.inv{D{i}.val-1};
            warning('Duplicating the last forward model for subject %d.', i);
        else
            error('Forward model is missing for subject %d.', i);
        end
    end
    
    D{i}.inv{D{i}.val}.inverse = inverse;
end

D = spm_eeg_invert(D);

if ~iscell(D)
    D = {D};
end

for i = 1:numel(D)
    save(D{i});
end

clear i inverse cfg job list mod P return_path sourcemodel2d transform;