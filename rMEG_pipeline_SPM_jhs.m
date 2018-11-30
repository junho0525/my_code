%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEG Pipeline using SPM
% 
% This Pipeline code is for SPM Example Dataset
% cdbespm12_SPM_CTF_MEG_example_faces1_3D.mat
%
% Last modified - 13/Dec/2017 by JHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
D = spm_eeg_load('cdbespm12_SPM_CTF_MEG_example_faces1_3D'); % load SPM_example preprocessed file
spm_eeg_inv_checkforward(D); % check forward model


%
% Simulation Code
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
job.isinversion.setsources.isSin.foi = [10, 20];
job.isinversion.setsources.dipmom=[10 10; 10 5];
job.isinversion.setsources.locs = [-6 28 34 ; -46 -70 36];
% [-6 28 34; -6 -52 40; -46 -70 36 ;-48 -69 7]; %ACC, PCC, LLP, LV1, LMT
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

clear dipfwhm dipmom fname fullname i j job meshsourceind mnimesh nAmdipmom Nsources ormni out simsignal sinfreq SNRdB timeind trialind whitenoisefT woi;

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


