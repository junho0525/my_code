%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This Is M/EEG Resting State Source Simulation & Localization Pipeline  %
%                           Code By Using SPM                             %
%                                                                         %
% HCP rest MEG data simulation & invert                                   %
%     100307_MEG_3-Restin_rmegpreproc                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally edited                                                          %
%     2017.12.13 xx:xx - By Junho Son                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load HCP Data and Preprocessing
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

%   High-pass Filter
S.D = D;
S.band = 'high';
S.freq = 0.1;
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
S.freq = 30;
D = spm_eeg_filter(S);
clear S;

%% Specify Forward Model
%opengl hardware;


return_path = pwd;

% ensure that the time and date of execution are not stored in the provenance information
%%%global ft_default
%%%ft_default.trackcallinfo = 'no';

load('100307_MEG_3-Restin_icaclass_vs.mat'); % This file is originally in Restin/icaclass directory of Restin_preproc

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
load('100307_MEG_anatomy_sourcemodel_2d.mat'); % This file is originally in anatomy directory
sourcemodel2d=ft_convert_units(sourcemodel2d, 'mm');
sourcemodelsubj = sourcemodel2d;

load('100307_MEG_anatomy_headmodel.mat');
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
D.inv = [];
D.inv{1}.forward = [];
D.inv{1}.method='Imaging';
D.inv{1}.forward.toMNI = transform.bti2spm;
%
D.inv{1}.forward.fromMNI = transform.spm2bti;
D.inv{1}.forward.voltype='Single Shell';
D.inv{1}.forward.vol = ft_convert_units(cfg.vol,'m');
D.inv{1}.forward.modality = 'MEG';
D.inv{1}.forward.siunits=1;
D.inv{1}.forward.mesh=[];
D.inv{1}.forward.mesh_correction=[];
D.inv{1}.forward.mesh.face = cfg.grid.tri;
grid = ft_convert_units(cfg.grid,'m');
D.inv{1}.forward.mesh.vert = grid.pos;
D.inv{1}.forward.sensors = [];
D.inv{1}.forward.sensors = ft_convert_units(cfg.grad,'m');
D.inv{1}.mesh.tess_mni.face = D.inv{1}.forward.mesh.face;
grid = ft_convert_units(cfg.grid,'mm');
D.inv{1}.mesh.tess_mni.vert = spm_eeg_inv_transform_points(D.inv{1}.forward.toMNI,grid.pos); % tess_mni should have mni coordinate
Datareg(1).modality = 'MEG';
D.inv{1}.datareg(1) = Datareg(1);
D.inv{1}.forward.toMNI = [D.inv{1}.forward.toMNI(:, 1:3)*1000 D.inv{1}.forward.toMNI(:,4)];
D.inv{1}.forward.fromMNI = eye(4)/D.inv{1}.forward.toMNI;
clear Datareg grid;
% check meshes and co-registration
spm_eeg_inv_checkforward(D); 

%% Simulation Code
%
% specify job for spm_cfg_eeg_inv_simulate.m, function run_simulation(job)
%
job = [];
job.D = [];
fname = strsplit(D.fname, '.');
job.D{1} = fname{1};
job.val = 1;
job.prefix = 'sim_';
job.whatconditions.all =1;
job.isinversion.setsources.woi = [100 400]; % Specify Time Window that source simulation starts and ends
% job.isinversion.setsources.isSin.foi = [10 ; 20]; % Specify Sinusoid Frequency
% job.isinversion.setsources.dipmom=[10 10; 10 5]; % Specify Dipole Moment
% job.isinversion.setsources.locs = [53.7203 -25.7363 9.3949; -52.8484 -25.7363 9.3949]; % Specify Source Location
job.isinversion.setsources.isSin.foi = [10 ; 20];
job.isinversion.setsources.dipmom=[10 10;10 5];
job.isinversion.setsources.locs = [53.7203 -25.7363 9.3949; -52.8484 -25.7363 9.3949];
% [-6 28 34; -6 -52 40; -46 -70 36 ;-48 -69 7; -39 -7 56; -28 52 19]; %ACC, PCC, LLP,
% LV1, LMT, LFEF, L aPFC
job.isSNR.setSNR = 0;

%
% Run Simulation
% 
% modified code of spm_cfg_eeg_inv_simulate.m
%
trialind=[];
if isfield(job.whatconditions, 'condlabel')
    trialind =D.indtrial( job.whatconditions.condlabel);
    if isempty(trialind),
        error('No such condition found');
    end;
end
if numel(job.D)>1,
    error('Simulation routine only meant to replace data for single subjects');
end;



if size(modality(D),1)>1,
    error('only suitable for single modality data at the moment');
end;

save(D);
D = {};


for i = 1:numel(job.D) %% only set up for one subject at the moment but leaving this for the future
    D{i} = spm_eeg_load(job.D{i});
    D{i}.val = job.val;
    
    D{i}.con = 1;
    if ~isfield(D{i}, 'inv')
        error(sprintf('Forward model is missing for subject %d', i));
    elseif  numel(D{i}.inv)<D{i}.val || ~isfield(D{i}.inv{D{i}.val}, 'forward')
        if D{i}.val>1 && isfield(D{i}.inv{D{i}.val-1}, 'forward')
            D{i}.inv{D{i}.val} = D{i}.inv{D{i}.val-1};
            warning(sprintf('Duplicating the last forward model for subject %d', i));
        else
            error(sprintf('Forward model is missing for subject %d', i));
        end
    end
    
end; % for i


if isfield(job.isSNR,'whitenoise'),
    whitenoisefT=job.isSNR.whitenoise; %% internal units as femto Tesla
    SNRdB=[];
else
    SNRdB=job.isSNR.setSNR;
    whitenoisefT=[];
end;

if isfield(job.isinversion,'setsources'), %% defining individual sources (This block should be executed)
    
    %%%
    Nsources=size(job.isinversion.setsources.locs,1)
    
    if (size(job.isinversion.setsources.dipmom,1)~=Nsources),
        error('Number of locations must equal number of moments specified');
    end;
    
    mnimesh=[]; %% using mesh defined in forward model at the moment
    
    
    ormni=[]; %% dipoles will get orientations from cortical mesh
    dipmom=job.isinversion.setsources.dipmom;
    if size(job.isinversion.setsources.dipmom,2)==3, %% dipole moment has 3 dimension
        disp('Simulating dipoles without reference to cortical surface');
        for j=1:size(dipmom,1),
            ormni(j,:)=dipmom(j,:)./sqrt(dot(dipmom(j,:),dipmom(j,:))); %% unit orientation in MNI space
            nAmdipmom(j)=sqrt(dot(dipmom(j,:),dipmom(j,:))); % magnitude of dipole
        end;
    else %% only one moment parameter given
        nAmdipmom=dipmom(:,1); %% total momnent in nAm
        dipfwhm=[];
        if size(dipmom,1)==2,
            dipfwhm=dipmom(:,2); %% fhwm in mm
        end;
    end;
    
    woi=job.isinversion.setsources.woi./1000;
    timeind = intersect(find(D{1}.time>woi(1)),find(D{1}.time<=woi(2)));
    simsignal = zeros(Nsources,length(timeind));
    
    if isfield(job.isinversion.setsources.isSin,'fromfile'),
        % Simulate orthogonal Gaussian signals
        filename=cell2mat(job.isinversion.setsources.isSin.fromfile);
        a=load(filename,'s');
        
        if size(a.s,1)~=Nsources,
            error('size of simulated data in file does not match number of sources');
        end;
        
        if size(a.s,2)~=length(timeind),
            warning(sprintf('sim signal from file is not same length as time window (%d vs %d samples) truncating or padding with zeros',size(a.s,2),length(timeind)));
        end;
        usesamples=1:min(length(timeind),size(a.s,2));
        simsignal(:,usesamples)=a.s(:,usesamples);
        
    end; % if isfield fband
    
    if isfield(job.isinversion.setsources.isSin,'fband'),
        % Simulate orthogonal Gaussian signals
        
        simsignal=randn(Nsources,length(timeind));
        % filter to bandwidth
        simsignal=ft_preproc_lowpassfilter(simsignal,D{1}.fsample,job.isinversion.setsources.isSin.fband(2),2);
        simsignal=ft_preproc_highpassfilter(simsignal,D{1}.fsample,job.isinversion.setsources.isSin.fband(1),2);
        [u,s,v]=svd(simsignal*simsignal');
        simsignal=u*simsignal; %% orthogonalise all signals
    end; % if isfield fband
    
    %  simsignal=simsignal.*repmat(nAmdipmom,1,size(simsignal,2)); %% now scale by moment
    
    if isfield(job.isinversion.setsources.isSin,'foi'),
        % simulate sinusoids
        sinfreq=job.isinversion.setsources.isSin.foi;
        % Create the waveform for each source
        
        for j=1:Nsources                % For each source
            simsignal(j,:)=sin((D{1}.time(timeind)- D{1}.time(min(timeind)))*sinfreq(j)*2*pi);
        end; % for j
        
        
    end; %% if isfield foi
    
    simsignal=simsignal./repmat(std(simsignal'),size(simsignal,2),1)'; %% Set sim signal to have unit variance
    %[D,meshsourceind,signal]=spm_eeg_simulate(D,job.prefix, job.isinversion.setsources.locs,simsignal,woi,whitenoisefT,SNRdB,trialind,mnimesh,SmthInit);
    figure();
    plot(D{1}.time(timeind),simsignal);
    xlabel('time');
    ylabel('Normalized amplitude (to be scaled by dip moment later)');
    legend(num2str([1:Nsources]'));
    [D,meshsourceind]=spm_eeg_simulate(D,job.prefix,job.isinversion.setsources.locs,simsignal,ormni,woi,whitenoisefT,SNRdB,trialind,mnimesh,dipfwhm,nAmdipmom);
    
else %% simulate sources based on inversion
    if ~isfield(D{i}.inv{job.val},'inverse'),
        error('There is no solution defined for these data at that index');
    end;
    
    [D]=spm_eeg_simulate_frominv(D,job.prefix,job.val,whitenoisefT,SNRdB,trialind);
end;



if ~iscell(D)
    D = {D};
end

for i = 1:numel(D)
    save(D{i});
end


fullname=[D{1}.path filesep D{1}.fname];
out.D = {fullname};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Invert SPM Example
%  spm_cfg_eeg_inv_invert.m -> function run_inversion(job)
%  Specify job
%

job = [];
job.D = D; % D is cell type
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

D = D{1};
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

